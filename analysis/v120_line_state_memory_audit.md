# SMElib line-state memory audit

## Scope

This audit investigates the high resident-memory peak of the v1.2.0 release
candidate on the canonical cool, metal-rich 800 A workload.  The run uses the
1,258,596-line input list, a 4250 K / log(g)=4.5 / [M/H]=+0.3 plane-parallel
atmosphere with 55 depth layers, `accrt=1e-4`, `accwi=3e-3`, ALMAX selection,
adaptive continuum opacity, and the production batched adaptive transfer grid.

Measurements are independent processes sampled at 10 ms.  The reported
"synthesis increment" is the peak reached after SME construction minus the RSS
immediately before synthesis.  It is more stable than subtracting absolute
process peaks because the Python allocator retains a variable amount of the
temporary VALD parser storage.

## Root cause

`InputModel` allocates three depth-major arrays for every input line:

```text
LINEOP[depth][line]
AVOIGT[depth][line]
VVOIGT[depth][line]
```

`ALMAXRange` calls `LINEOPAC` for every line and commits those pages.  For this
workload, the double-precision logical payload is:

```text
1,258,596 lines * 55 depths * 3 arrays * 8 bytes
= 1,661,346,720 bytes (1.547 GiB)
```

The observed RSS jump at `ALMAXRange` agrees with that payload.  The adaptive
transfer output is not the dominant allocation: its reserved capacity is about
259 MiB and its 532,481 actual nodes contain about 61 MiB of wave/intensity
data.

The benchmark harness's historical deep copy of the pandas line list is a
secondary cost.  It changes the measured peak from about 4.86 GB to 5.37 GB,
but it does not explain the native `ALMAXRange` jump.

There is a separate parser peak.  `ValdFile.loads` currently holds the complete
file as Python strings, creates slices, joins the line records into another
large string, and then asks pandas to parse that string.  The parser reached
about 4.19 GB RSS in the mixed-precision run.  Consequently, eliminating all
native line-state storage would still leave the end-to-end process peak near
the parser peak until the VALD reader is made streaming.

## Mixed-precision prototype

The prototype stores `LINEOP`, `AVOIGT`, and `VVOIGT` as `float`, while all
opacity accumulation, Voigt evaluation, and transfer arithmetic remains
double precision.  This changes cache storage precision, not the physical
algorithm or selection threshold.

| Metric | double baseline | float cache | Difference |
|---|---:|---:|---:|
| Logical three-array payload | 1,661,346,720 B | 830,673,360 B | -830,673,360 B |
| Post-build to synthesis peak | about 1.90 GB | about 1.09 GB | about -0.81 GB |
| Sampled absolute peak | 4,860,559,360 B | 4,258,267,136 B | -602,292,224 B |
| Synthesis wall time | 185.6 s | 180.4 s | no measured slowdown |

The smaller change in absolute peak is caused by different retained parser
memory in the two independent processes.  The stage-local reduction matches
the theoretical cache saving.

### Numerical comparison

The 800 A spectrum was compared point-for-point with the committed release
candidate double-cache result on the identical 0.05 A output grid:

```text
max |delta normalized flux| = 1.3121e-8
RMS delta flux              = 1.2325e-9
99th percentile |delta F|   = 4.5457e-9
delta EW over 800 A         = -6.8681e-8 A
selection-mask changes      = 0
physical-range changes      = 0
```

Ten-Angstrom solar, cool metal-rich dwarf, and moderate giant cases had maxima
of `6.6e-9` to `9.3e-9`.  Additional metal-poor, 3500 K spherical giant, Mg b,
and H-alpha cases had maxima of `2.5e-9` to `1.6e-8`, again with no mask or
range changes.

The complete source test suite passes:

```text
187 passed, 1 xfailed, 15 warnings
```

## Do all three arrays need to remain resident?

Not equally.

`VVOIGT` is algebraically redundant.  Its non-hydrogen value is the inverse
Doppler width,

```text
1 / (lambda_0 * sqrt(T_depth * mass_factor_species + vturb_depth^2)),
```

and its hydrogen value is the corresponding dimensionless Doppler width.  It
can therefore be represented by a small `depth * species` table plus one
inverse central wavelength per line.  Removing the full float `VVOIGT` matrix
would save about 264 MiB on this workload without introducing a physical
approximation.  This is the preferred next low-risk prototype.

`LINEOP` and `AVOIGT` are different.  The current wavelength-centric `OPMTRX`
reads them for every active line at every transfer node and depth.  Replacing
those reads with direct on-demand calls would repeatedly recompute populations,
exponentials, damping terms, and profile parameters, and is expected to cause a
large runtime regression.

The promising structural alternative is wavelength-chunked materialization:
retain immutable ALMAX masks and physical ranges globally, but construct line
state only for lines whose support overlaps the current wavelength chunk.  In
the cool 800 A result, a 50 A chunk overlaps at most 74,557 selected lines.  Its
three float matrices would be about 47 MiB, versus 792 MiB for all input lines,
and a line overlaps only 1.012 chunks on average.  This makes chunking worth a
bounded prototype, but it requires changes to adaptive-node scheduling and
should not be folded into the mixed-precision patch.

## Recommendation

1. Keep the all-float cache change as the small, independently testable memory
   fix, subject to normal review.
2. Prototype factorized `VVOIGT` separately and require no meaningful runtime
   regression plus the same numerical checks.
3. Treat chunked `LINEOP`/`AVOIGT` materialization as a larger performance task
   with a go/no-go benchmark, not as release-candidate cleanup.
4. Audit the VALD parser separately; it becomes the absolute RSS ceiling once
   native synthesis storage is reduced.
