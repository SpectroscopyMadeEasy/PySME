# PySME v1.1.1 to current: implementation and runtime summary

## Scope

This note summarizes the current `codex/cdr-physical-range-selection` branch at
`394e587` relative to PySME `v1.1.1` (`56939dd`).  The bundled SMElib revisions
are `d80f524` and `533964f`, respectively.

The comparison covers both user-visible default behavior and component-level
benchmarks.  It is important not to present a continuum-cache speedup, an
abundance-fit state-reuse speedup, and their combined end-to-end speedup as if
they were three independent multiplicative results.

## Code changes

### Exact EOS history reuse

- Abundance fitting enables a scoped SMElib EOS warm start.
- The previous converged electron pressure and species partial pressures are
  used only as the next iteration's initial guess.
- Model layers, species, temperature, pressure, abundance scaling, EOS mode,
  and finite/positive values are validated before reuse.
- The full EOS solve still runs, so this is an iteration-count optimization,
  not a physical approximation.

### Serialized SMElib workflows

- Python entry points that touch SMElib's process-global native state are
  serialized with a reentrant session lock.
- This fixes races between multiple synthesizer/solver objects in one process.
- It deliberately prevents unsafe same-process threaded overlap; process-based
  line-info workers remain available.

### Native ALMAX and line-validity pipeline

- The default line-selection method changed from `internal` to `almax`.
- SMElib can receive and reuse precomputed line ranges, strong-line masks, and
  central depths.
- Fixed-grid transfer now discovers opacity-based validity ranges instead of
  returning the `wlcent +/- 150 A` placeholders installed by `InputLineList`.
- An interval sweep restricts wavelength-point candidate lookup to overlapping
  validity ranges, with an automatic cost check that keeps the full scan when
  the index would not save enough work.
- The immediately following first `Transf` can consume line opacity and Voigt
  state computed by `ALMAXRange`, avoiding a duplicate `LINEOPAC` pass.
- Cumulative wavelength-bin strong-line selection moved into SMElib so CDR and
  optional binned-ALMAX workflows use the same native implementation.

### Adaptive continuum-opacity grid

- The default continuum mode changed from exact `CONTOP` at every query to an
  adaptive grid with nominal 1 A spacing.
- Known H I thresholds and Mg I/Si I PEACH knots are inserted explicitly.
- A `+/-0.02 A` exact guard band protects sharp edge shoulders; other intervals
  are recursively refined from curvature probes down to a configurable minimum
  step.
- All 13 continuum source components are cached, preserving consistent true
  absorption, coherent scattering, total extinction, and scattering source.
- Public controls are `continuum_grid`, `continuum_grid_base_step`,
  `continuum_grid_rtol`, and `continuum_grid_min_step`; `"exact"` remains a
  continuum-only diagnostic/reference mode.
- A redundant outer `CentralDepth -> CONTOP` call was removed.

### Prepared abundance synthesis

- The first abundance-only `solve(...)` evaluation installs the complete line
  list and atmosphere/model in SMElib.
- Later evaluations reuse those inputs, but an abundance change invalidates
  ALMAX because line strength and validity range are abundance-dependent.  The
  repeated sequence is therefore
  `InputAbund -> Ionization -> Opacity -> ALMAXRange ->`
  `InputLinePrecomputedInfo -> Transf`.
- ALMAX is evaluated after the new abundance reaches SMElib.  Its immediately
  following `Transf` reuses the resulting line-opacity/Voigt state, so the
  correctness fix does not add a second `LINEOPAC` pass.
- The effective 99-element abundance vector is recorded with ALMAX metadata;
  any element change invalidates it.  A stale atmosphere-only ALMAX cache is
  bypassed after an observed abundance change.
- Prepared abundance solving is limited to native `almax` and `internal`
  selection.  CDR selection falls back to full setup because its metadata is
  not yet abundance-aware.
- Mixed physical-parameter fits, dynamic line-list mode, forced line-info
  recomputation, callable atmospheres, and line-list parameter fits retain the
  full setup path.
- The prepared path is bit-for-bit identical to the current full-reinput path.

### Correctness and maintenance changes

- NumPy record scalars from atmosphere data are handled correctly.
- Native/Python wrappers, C-wrapper entry points, changelog, reference docs,
  and regression tests were extended for the new state and selection APIs.
- The explored line-centric opacity-loop rewrite was removed after it failed
  the end-to-end stop condition; production retains one opacity engine.

Relative to v1.1.1, the parent repository changes 29 files by about 3,475
insertions and 140 deletions.  The SMElib submodule changes nine files by about
1,709 insertions and 70 deletions.  Current verification is 174 passed and one
expected failure.

## Version-to-version fixed-grid runtime

### Benchmark definition

- Machine: Apple M5, 10 cores, macOS 26.6.2 arm64, Python 3.11.14.
- Separate isolated wheels were built from v1.1.1 and current `394e587`.
- Line list: `merged_3700-9500_hfs.lin`, extracted to a 10 A output window at
  5200 A plus about 150 A of line-list padding: 378,154 input lines.
- Fixed transfer/output grid: `sme.wave = sme.wint`, 0.05 A, 201 points.
- `accrt=1e-4`; VALD parsing is excluded.
- Each version uses its own defaults.  Current first-call time includes ALMAX
  and physical-range construction.  The repeated current call uses the
  automatic abundance-only prepared state.  The repeated v1.1.1 call performs
  its normal full reinput.
- The two-call workload is one initial synthesis at the baseline Fe abundance
  followed by one `+0.05 dex` Fe evaluation.

| model | v1.1.1 initial | current initial | initial speedup | v1.1.1 repeat | current repeat | repeat speedup | two-call speedup |
|---|---:|---:|---:|---:|---:|---:|---:|
| solar dwarf | 47.073 s | 2.056 s | 22.89x | 47.472 s | 1.915 s | 24.79x | 23.81x |
| metal-poor dwarf | 44.427 s | 1.990 s | 22.33x | 45.047 s | 1.881 s | 23.95x | 23.12x |
| cool metal-rich dwarf | 53.602 s | 3.006 s | 17.83x | 52.952 s | 2.729 s | 19.40x | 18.58x |

The repeat speedup is the combined result of adaptive continuum evaluation,
ALMAX/range selection, EOS warm-start support, and prepared state reuse.  It is
not the isolated effect of prepared synthesis.  Unlike the earlier prototype,
these measurements recompute ALMAX whenever Fe changes.  Prepared state alone
is consequently a small setup optimization, not the source of a 100x repeat
speedup.

For a 200 A solar window with 592,757 input lines and 4,001 output points,
current completed the initial default synthesis in 3.442 s.  The v1.1.1 run
was stopped after 360 s, giving a conservative lower-bound speedup of more
than 104x for this deliberately large full-line-list workload.

These ratios are workload-dependent.  They apply to the large HFS line list
with wide line padding requested for this investigation.  A user who supplies
an already narrow, aggressively filtered line list will see a smaller ratio.

## Version-to-version wave-only runtime

The normal `sme.wave`-only case was measured separately.  In this mode
`sme.wint` is not set and direct synthesis uses
`reuse_wavelength_grid=False`, so every call lets SMElib construct its own
adaptive transfer grid before PySME interpolates onto the 201-point output
grid.  All other inputs match the fixed-grid comparison.

| model | v1.1.1 initial | current initial | initial speedup | v1.1.1 repeat | current repeat | repeat speedup | two-call speedup |
|---|---:|---:|---:|---:|---:|---:|---:|
| solar dwarf | 34.794 s | 3.371 s | 10.32x | 35.875 s | 3.222 s | 11.14x | 10.72x |
| metal-poor dwarf | 33.258 s | 2.286 s | 14.55x | 33.164 s | 2.045 s | 16.22x | 15.34x |
| cool metal-rich dwarf | 230.196 s | 116.500 s | 1.98x | 225.559 s | 116.564 s | 1.94x | 1.96x |

The initial default spectra agree much more closely in this adaptive-transfer
case than in the fixed-grid comparison: maximum absolute normalized-flux
differences are `7.82e-10`, `7.16e-10`, and `1.40e-9` for the solar,
metal-poor, and cool models.  This confirms that the larger fixed-grid version
difference comes mainly from the changed fixed-grid line-range/weak-line
semantics.

After changing Fe by `+0.05 dex`, maximum version differences are
`1.72e-7`, `5.60e-8`, and `6.11e-8`.  Current now recomputes ALMAX/ranges for
the changed abundance instead of retaining the line set from the first call.
This removes the previous `1e-4`-level repeat discrepancy while keeping the
current wave-only path about 2--16x faster for these models.

The wave-only result materially changes the performance interpretation.  Solar
and metal-poor synthesis still improve by roughly 10--16x, but the line-rich
cool model improves by only about 1.9x.  The fixed-grid interval sweep is not
used in this mode, and the cool adaptive transfer grid itself remains expensive.

## End-to-end fixed-grid `solve()` runtime

The version comparison also measures a real one-parameter abundance fit, not
only direct calls to the synthesizer:

- Each revision fits its own noise-free default spectrum.  This isolates
  optimizer runtime from the already documented default-spectrum difference
  between revisions.
- The fit starts with Fe displaced by `+0.10 dex` and uses an uncertainty of
  `0.01` at every wavelength point.  All least-squares stopping criteria retain
  their normal defaults.
- For the current ALMAX path, target construction exercises the same abundance
  transition as the fit and recomputes ALMAX at the target abundance.
- Target construction and VALD parsing are excluded.  The wavelength window,
  378,154-line input, grid, atmosphere models, and `accrt=1e-4` are the same as
  in the 10 A synthesis benchmark above.

| model | v1.1.1 solve | current solve | residual iterations (old/new) | synthesis calls (old/new) | time per synthesis (old/new) | per-call speedup | solve speedup |
|---|---:|---:|---:|---:|---:|---:|---:|
| solar dwarf | 493.10 s | 18.48 s | 5 / 5 | 10 / 10 | 49.31 / 1.846 s | 26.71x | 26.69x |
| metal-poor dwarf | 442.49 s | 18.34 s | 5 / 5 | 10 / 10 | 44.25 / 1.832 s | 24.15x | 24.13x |
| cool metal-rich dwarf | 323.64 s | 22.66 s | 3 / 4 | 6 / 8 | 53.94 / 2.831 s | 19.05x | 14.28x |

The solar and metal-poor comparisons have identical optimizer and synthesis
call counts, so their approximately 24--27x end-to-end ratios directly reflect
the reduced cost per evaluation.  The cool model takes two additional synthesis
calls in the current revision; its 14.28x total speedup is consequently smaller
than its 19.05x per-call speedup.  All three current fits recover the target Fe
abundance to within `1.6e-7 dex`; the largest final absolute flux residual is
`1.82e-7`.

This is the combined default-path benefit of the line-selection/range work,
adaptive continuum evaluation, EOS history, and prepared abundance state.  It
must not be quoted as the isolated contribution of Prepared Synthesis.  A
200 A v1.1.1 solve was not attempted because even its first synthesis exceeded
360 s in the direct benchmark; a multi-evaluation solve would not add a useful
constraint beyond that existing lower bound.

## End-to-end wave-only `solve()` runtime

The abundance solve was also rerun without setting `sme.wint`.  Target spectra
are generated with the same revision and transfer mode.  Target construction
starts at `+0.10 dex` and then recomputes ALMAX at the target abundance, matching
the lifecycle used during the fit.

| model | current wave-only solve | current fixed-grid solve | wave-only / fixed cost | residual iterations | synthesis calls |
|---|---:|---:|---:|---:|---:|
| solar dwarf | 25.74 s | 18.48 s | 1.39x | 5 | 10 |

The post-fix solar fit recovers the target Fe abundance within `1.3e-7 dex`.
Metal-poor and cool wave-only solves were not repeated after the ALMAX
invalidation correction; their earlier values used a fixed ALMAX line set and
must not be quoted as current results.

The v1.1.1 solar wave-only solve was stopped after 1,011 s without completing,
at the user's requested stop condition.  Against the completed current time of
25.74 s, this establishes a conservative end-to-end speedup lower bound of
more than 39x.  The v1.1.1 metal-poor and cool wave-only solves were not run
after that stop condition was reached.  The incomplete solar run produced no
final iteration/call count, so no per-synthesis normalization is claimed for
this lower bound.

## Component-level runtime

- Production adaptive continuum grid: the directly comparable solar
  `ALMAXRange + CentralDepth` chain fell from 61.23 s to 3.56 s, or 17.2x.
- The production grid reduces roughly 0.76--1.14 million continuum queries to
  about 1,246--1,420 exact continuum nodes in the tested CDR workload.
- Prepared state still avoids repeated line-list/model input, but ALMAX must be
  recomputed after abundance changes; its isolated saving is therefore small
  compared with the earlier fixed-ALMAX prototype.
- A numerically equivalent line-centric loop transpose improved complete cool
  synthesis by only 1.12x at 200 A and 1.18x at 800 A, so it was not retained.

## Numerical behavior relative to v1.1.1

Default-to-default normalized-flux differences for the 10 A fixed-grid version
benchmark were:

| model | maximum absolute difference | RMS difference |
|---|---:|---:|
| solar dwarf | `3.59e-4` | `2.21e-4` |
| metal-poor dwarf | `1.17e-4` | `5.32e-5` |
| cool metal-rich dwarf | `1.10e-3` | `5.05e-4` |

These differences must not be attributed to continuum interpolation alone.
The adaptive-continuum production benchmark found a worst tested normalized
flux difference of `7.24e-6` in the cool cumulative-selection boundary case,
with most models near `1e-9`; the Mg I 3756.6 A edge regression is
`3.58e-11` after adding the exact guard band.

The larger version-to-version difference comes mainly from changed line-range
and weak-line selection semantics.  v1.1.1 fixed-grid synthesis retained many
weak contributions inside broad placeholder ranges.  Current synthesis uses
opacity-based ranges and default ALMAX selection at the local `accrt` ratio.
Many individually weak lines can accumulate into a visible pseudo-continuum,
so `accrt=1e-4` is not a `1e-4` bound on final normalized flux.

Setting `continuum_grid="exact"` only restores exact continuum evaluation; it
does not restore v1.1.1 fixed-grid line-range behavior.  In the solar 10 A
test, current `internal + exact continuum` still differed from v1.1.1 by about
`3.58e-4`, and took about 60 s versus 47 s.  It should therefore be documented
as a continuum reference mode, not as a complete v1.1.1 compatibility mode.

## Documentation plan

1. **Release notes / changelog**: summarize the changed defaults, headline
   version-to-version runtime, prepared abundance fitting, thread-safety fix,
   and the numerical-compatibility caveat.
2. **Line-selection reference**: keep the detailed ALMAX/CDR/range semantics,
   `accrt` definition, interval index, and continuum-grid controls here.
3. **Performance and accuracy page**: add a stable user-facing table with the
   benchmark inputs, hardware, initial versus repeated timing, and flux
   differences.  Link to the analysis drivers rather than embedding every
   development experiment.
4. **Migration note from v1.1.1**: explicitly state that default spectra are
   not bit-identical, `continuum_grid="exact"` is not a full legacy mode, and
   users requiring frozen legacy results should retain v1.1.1 until a complete
   legacy line-range option exists or revalidate their workflow.
5. **Developer analysis**: retain continuum, prepared-state, and rejected
   line-centric reports under `analysis/`; do not move their raw NPZ/JSON data
   into normal user documentation.
