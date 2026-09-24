# Spherical generation-batched RKINTS audit and benchmark

## 1. Problem and conclusion

The spherical adaptive transfer path had the same wavelength-scheduling
pathology as plane-parallel RKINTS: it evaluated one wavelength at a time, so
`OPMTRX` could not use a sorted-generation interval sweep.  The spherical
formal solution itself was not the bottleneck.

The production candidate now uses the same generation scheduler for
plane-parallel and spherical atmospheres.  Geometry-specific evaluators are
limited to the final formal solution:

```text
shared immutable ALMAX mask + physical Wlim
        -> shared adaptive generation scheduler
        -> LineIntervalSweep + OPMTRXIndexed
        -> plane-parallel: TBINTG
        -> spherical: ray path + grazing mirror + TBINTG_sph
```

The implementation is a **strong GO**.  Controlled spherical intensities are
bitwise identical to fixed-grid evaluation at every accepted node and ray.
The three requested dense-window MARCS regimes accelerate by 31.35x, 23.88x,
and 6.88x; Mg b and H-alpha accelerate by 39.46x and 49.83x.

## 2. Source audit of the legacy spherical path

The relevant source is `smelib/src/sme/sme_synth_faster.cpp`.

- `Transf` starts at line 7369.  It validates line state, selects the largest
  `mu` as `IMU_REF`, restores or computes line opacity/ranges, and builds the
  active-line index.
- Spherical geometry is constructed from line 7689.  For each ray it computes
  the impact parameter, a ray-specific column-mass path, `NRHOXs[imu]`, and a
  `grazing[imu]` flag.
- `BuildSphericalRayOrder` is at line 8162 and `RKINTS_sph` starts at line
  8174.
- The shared generation scheduler starts at line 8625.  The spherical batch
  evaluator starts at line 8567 and its thin `RKINTS_sph_batched` wrapper at
  line 8765.
- `TBINTG_sph` starts at line 9257.
- `OPMTRXIndexed` and `OPMTRX` start at lines 10478 and 10489.

Legacy spherical initial-grid construction is the same basic geometry as
plane-parallel RKINTS:

1. Evaluate the blue endpoint.
2. For every active line centre separated from the previous node by more than
   0.3 km/s, insert the preceding midpoint and the line centre.
3. Evaluate the red endpoint.
4. Walk left to right, insert every tested midpoint, and recursively refine
   intervals that fail the intensity criterion.

The criterion is evaluated on the largest-`mu` (disk-centre) ray:

```text
error = (|I_mid - (I_left + I_right)/2|
         + 0.005 |I_left - I_right|) / I_cont(mid)
```

An interval is accepted when `error < accwi`, or when its left half-width is
at the 0.3 km/s minimum.  Spherical legacy RKINTS uses the left endpoint to
convert that velocity floor to wavelength; the shared scheduler preserves
this small difference from the current plane-parallel batched convention.

Unlike legacy plane-parallel RKINTS, legacy `RKINTS_sph` does **not** apply a
line-centre `MARK=2` second-pruning step.  It did, however, shorten
`Wlim_left/right` during its midpoint/line-centre seed scan when the sampled
`ALMAX` was below `accrt`.  The new spherical path changes neither the active
mask nor the physical ranges.

## 3. Geometry-specific radiative transfer

The wavelength-independent spherical geometry is built once per `Transf`:

- Every `mu` has a different impact parameter and column-mass path.
- Normal rays traverse all `NRHOX` layers with a depth-dependent geometric
  path factor.
- Grazing rays stop at their deepest intersected shell.  Their path is then
  mirrored across the tangent point, so `NRHOXs[imu]` is twice the number of
  near-side shells reached.
- The radial opacity and source arrays from `OPMTRX` are mirrored into the
  far-side half before `TBINTG_sph`.
- `TBINTG_sph` uses a zero incoming boundary intensity for a grazing ray and
  the usual deep boundary approximation for a normal ray.

The quantities depending on wavelength are the opacity/source arrays and the
emergent intensities.  The ray paths, layer counts, grazing flags, and
far-side mapping do not depend on wavelength and are reused for every batch
probe.

Therefore generation batching changes only wavelength scheduling.  It does
not change the spherical formal-transfer equation.

### Grazing-ray order correction

The audit exposed one pre-existing implementation issue.  Grazing rays mirror
values in-place in shared work arrays.  If a grazing ray was evaluated before
a normal ray, it could overwrite depth entries subsequently consumed by that
normal ray, making intensities depend on caller `mu` order.

The correction is deliberately minimal: for each wavelength, normal rays are
evaluated first, followed by grazing rays from deepest to shallowest.  Results
are still written into the original caller ray slots.  This prevents an
in-place far-side write from contaminating any later near-side input, without
changing `TBINTG_sph` or adding a second opacity engine.  The same local ray
ordering is used by legacy/fixed spherical transfer, batched transfer, and the
spherical contribution-function path.

Forward and reversed `mu` arrays now give bitwise-identical wavelengths,
line intensities, and continuum intensities after restoring ray order.

## 4. Reuse of the indexed opacity path

`LineIntervalSweep` depends only on `MARK`, `Wlim_left/right`, and a sorted
wavelength array.  It has no atmosphere-geometry dependency.  Spherical
fixed-grid synthesis already used `LineIntervalSweep + OPMTRXIndexed`, which
provided a direct proof that no `OPMTRXIndexed_sph` was required.

`OPMTRX` constructs line/continuum opacity and source arrays on the radial
atmospheric depth scale.  All rays consume those same radial arrays; only the
subsequent path mapping and formal solution are ray-specific.  Indexed
candidates remain sorted in original line-list order, so opacity accumulation
order is unchanged.

Hydrogen profiles, autoionization checks, ordinary Voigt profiles, continuum
scattering, and NLTE departure-coefficient opacity/source terms all remain
inside the existing `OPMTRX` implementation.  Batching does not bypass them.

## 5. Shared implementation

`RunBatchedAdaptiveRKINTS` owns all shared scheduling logic:

- immutable line-centre seed construction;
- complete generation-0 evaluation;
- sorted midpoint generations;
- endpoint/result caching with no recomputation;
- the `accwi` emergent-intensity test;
- final wavelength sorting and output layout.

It accepts a synthesis-local evaluator context and function pointer.  The
plane-parallel evaluator calls `OPMTRXIndexed + TBINTG`; the spherical
evaluator calls `OPMTRXIndexed`, applies the existing ray/far-side mapping,
and calls `TBINTG_sph` for line and continuum.

No new process-global wavelength or ray cache was introduced.  Node arrays,
generation arrays, interval lists, ray geometry pointers, and result arrays
are local to the current `Transf`.  Execution continues to use PySME's
existing process-wide SMElib lock.

With `SME_TIMING=1`, the scheduler reports aggregate generation probes and
accepted/refined intervals in addition to the existing opacity, formal-RT,
and interval-candidate counters.

## 6. Production gating

Spherical batching is selected only when all of the following hold:

```text
MOTYPE == spherical
NWL == 0                         # adaptive grid, not a supplied fixed grid
long_continuum == true
transfer_grid_method == batched
validated precomputed mask/ranges are active
```

The following still use their established paths:

- a user-supplied fixed wavelength grid;
- internal line selection or missing/stale precomputed line information;
- `long_continuum=False`;
- explicit `transfer_grid_method="legacy"`.

Internal line selection was not changed by this task.

## 7. Benchmark setup

All main cases use real spherical MARCS model files with 56 depth layers:

| model | Teff | log g | [M/H] | MARCS geometry |
|---|---:|---:|---:|---|
| moderate giant | 4500 K | 2.0 | 0.0 | spherical |
| cool line-rich giant | 3500 K | 1.0 | 0.0 | spherical |
| metal-poor giant | 4500 K | 2.0 | -2.0 | spherical |

Other settings:

- `accrt = 1e-4` for the immutable ALMAX membership mask;
- physical ranges computed with `line_select_range_floor = 1e-6`;
- `accwi = 3e-3`;
- seven Gauss-Legendre `mu` values from 0.02545 to 0.97455;
- the outer two rays are grazing in all three MARCS models;
- `long_continuum=True`;
- single-process wall times on the same machine;
- timings exclude VALD parsing and ALMAX/range discovery, but include the
  normal `Transf` line-opacity setup.

The dense and Mg b cases use 378,154 lines extracted with about 150 A of
buffer.  The H-alpha case uses 605,873 buffered lines.  The dedicated NLTE
case uses the tested 983-line H-alpha regression list and the cached
`nlte_H_pysme.grd` departure grid.

## 8. Grid and native-intensity equivalence

For all five LTE cases, historical legacy and new batched scheduling happened
to produce the exact same wavelength nodes:

| case | legacy nodes | batch nodes | common | legacy-only | batch-only | Hausdorff |
|---|---:|---:|---:|---:|---:|---:|
| dense moderate | 6,553 | 6,553 | 6,553 | 0 | 0 | 0 A |
| dense cool | 6,665 | 6,665 | 6,665 | 0 | 0 | 0 A |
| dense metal-poor | 863 | 863 | 863 | 0 | 0 | 0 A |
| Mg b moderate | 16,649 | 16,649 | 16,649 | 0 | 0 | 0 A |
| H-alpha moderate | 13,497 | 13,497 | 13,497 | 0 | 0 | 0 A |

The controlled reference has stronger meaning than the historical comparison:

- an immutable-state depth-first sequential replay produced exactly the same
  node array as generation batching;
- evaluating all accepted batch nodes through fixed-grid spherical `Transf`
  produced bitwise-identical `wave`, `sint`, and `cint`;
- equality holds for every one of the seven rays, including the two grazing
  rays, intermediate rays, and the disk-centre reference ray;
- the batch changed zero active physical-range values in every case.

Legacy shortened 2, 31, 4, and 4 active range endpoints in the dense
moderate, dense metal-poor, Mg b, and H-alpha cases respectively.  It changed
none in the dense cool case.  These historical range mutations explain the
small legacy-versus-new intensity differences below; batching itself is exact
under controlled semantics.

## 9. Flux and equivalent-width validation

Native intensities were first resampled ray-by-ray to a shared regular
log-wavelength grid at 0.3 km/s.  Only then were projected-area integration
and instrumental convolution applied.  No irregular grid was passed into the
index-based broadening/integration routine.

| case | max intrinsic delta | intrinsic RMS | max delta R=20k | max delta R=60k | EW delta |
|---|---:|---:|---:|---:|---:|
| dense moderate | 8.61e-8 | 4.56e-9 | 8.17e-9 | 2.24e-8 | 3.78e-9 A |
| dense cool | 0 | 0 | 0 | 0 | 0 A |
| dense metal-poor | 6.43e-6 | 3.55e-7 | 5.09e-7 | 1.20e-6 | 5.67e-7 A |
| Mg b moderate | 6.52e-8 | 2.06e-9 | 4.88e-9 | 1.37e-8 | 3.35e-9 A |
| H-alpha moderate | 9.18e-8 | 3.65e-9 | 7.73e-9 | 2.16e-8 | 7.56e-9 A |

These are historical-legacy versus new immutable-range semantics, not errors
from generation batching.  The controlled fixed-node differences are exactly
zero.

## 10. Performance and bottleneck

| case | nodes | legacy Transf | batched Transf | speed-up |
|---|---:|---:|---:|---:|
| dense moderate | 6,553 | 108.336 s | 3.456 s | **31.35x** |
| dense cool | 6,665 | 132.959 s | 5.567 s | **23.88x** |
| dense metal-poor | 863 | 0.355 s | 0.052 s | **6.88x** |
| Mg b moderate | 16,649 | 280.284 s | 7.103 s | **39.46x** |
| H-alpha moderate | 13,497 | 374.045 s | 7.507 s | **49.83x** |

The number of physical wavelength evaluations is unchanged.  The benefit is
candidate reduction:

| case | legacy OPMTRX | batch OPMTRX | batch candidate/full visits | reduction | batch TBINTG_sph |
|---|---:|---:|---:|---:|---:|
| dense moderate | 107.154 s | 2.328 s | 5.89 M / 2.172 B | 99.73% | 0.040 s |
| dense cool | 131.749 s | 4.201 s | 10.79 M / 2.518 B | 99.57% | 0.037 s |
| dense metal-poor | 0.310 s | 0.016 s | 0.037 M / 6.60 M | 99.44% | 0.005 s |
| Mg b moderate | 279.017 s | 5.824 s | 14.73 M / 5.519 B | 99.73% | 0.103 s |
| H-alpha moderate | 372.002 s | 5.131 s | 11.03 M / 8.108 B | 99.86% | 0.083 s |

`TBINTG_sph` remains far below one percent of legacy wall time in the rich
cases.  After batching, `OPMTRX` is still the largest transfer component;
spherical ray integration does not become the dominant bottleneck.

The initial grid plus one midpoint generation was sufficient in all rich-line
cases.  The metal-poor case used six generations, 863 total unique probes, 104
refined intervals, and 431 accepted intervals.  Existing endpoints were not
recomputed.

## 11. Comparison with plane-parallel batching

The mechanism is the same in both geometries: nearly the same number of
`OPMTRX` calls is retained, while interval indexing removes almost all
irrelevant line visits.

| representative 10 A case | plane-parallel speed-up | spherical speed-up |
|---|---:|---:|
| solar/moderate metallicity | 10.24x solar dwarf | 31.35x moderate giant |
| line-rich cool | 23.91x cool dwarf | 23.88x cool giant |
| metal-poor | 5.85x dwarf | 6.88x giant |
| giant | 23.53x plane-parallel giant | 31.35x spherical giant |

The larger spherical moderate result is not a new ray optimization.  Its
legacy `OPMTRX` full-scan cost is especially severe, while its additional
`TBINTG_sph` work remains small.

## 12. NLTE, H lines, and special profiles

Architecture audit shows that batch probes enter the same `OPMTRX` code as
legacy/fixed synthesis:

- NLTE extinction and source terms use `BNLTE_low/BNLTE_upp` inside
  `OPMTRX`;
- H lines call the existing `hlinprof_` branch;
- autoionization validity checks and ordinary Voigt accumulation are
  unchanged;
- candidate filtering still uses the precomputed physical Wlim intervals.

Two NLTE validations passed:

1. A regression test installs non-unity departure coefficients in a
   spherical atmosphere with grazing rays.  Batched adaptive and fixed-node
   `sint/cint` are bitwise identical.
2. A real spherical 4500/2.0 MARCS atmosphere with the cached tested
   `nlte_H_pysme.grd` H-alpha departure grid produced 1,633 identical legacy
   and batch nodes.  Batched and fixed-node line/continuum intensities were
   bitwise identical, as were historical legacy and batch fluxes.  The small
   983-line case accelerated only 1.35x (0.198 s to 0.147 s), as expected when
   interval pruning removes much less work.

The buffered LTE H-alpha test separately validates the expensive H-profile
path at scale and accelerated by 49.83x.  Autoionization was audited in source
but not isolated as a dedicated empirical benchmark; it remains a documented
coverage limitation rather than an architectural bypass.

## 13. Regression tests

The regression suite now covers:

- spherical and plane-parallel immutable mask/range behavior;
- first-call and `keep_lineop` adaptive reuse;
- batched adaptive versus exact fixed-node line/continuum intensities;
- grazing, intermediate, and disk-centre rays;
- explicit spherical legacy fallback and historical range mutation;
- spherical non-LTE departure coefficients;
- interval-index exactness on fixed grids;
- forward/reversed `mu` order equality for spherical intensities;
- existing continuum-scattering and high-level convolution/integration tests.

Full result:

```text
185 passed, 1 xfailed, 15 warnings in 48.22 s
```

No brittle wall-time assertion was added to CI.  Path use is proven by the
immutable-range and fixed-node equality tests; performance is recorded in the
benchmark artifacts.

## 14. Known limitations

- Production batching still deliberately excludes internal line selection,
  fixed grids, missing/stale precomputed line state, and
  `long_continuum=False`.
- `accwi` remains the historical disk-centre emergent-intensity midpoint
  heuristic, not a global final-flux error bound.
- The large benchmark timings are single runs because each historical rich
  case costs minutes.  The native phase decomposition and consistent
  candidate-reduction ratios make the mechanism unambiguous, but the last
  digits of wall-time speed-ups should not be treated as statistically exact.
- A dedicated autoionization-only spherical spectrum was not run.

## 15. Final answers

1. **How similar are the schedulers?**  Seed geometry and intensity refinement
   are essentially the same.  Spherical legacy lacked plane-parallel's
   `MARK=2` second pruning, but did mutate Wlim; its minimum-spacing wavelength
   uses the left endpoint.
2. **Can the scheduler be shared?**  Yes.  One shared scheduler now serves
   both geometries through small evaluator contexts.
3. **Do `LineIntervalSweep/OPMTRXIndexed` apply to spherical transfer?**  Yes;
   line candidates and radial opacity/source construction are geometry
   independent.  No spherical opacity engine is needed.
4. **Are native ray intensities equivalent?**  Yes, bitwise under controlled
   immutable semantics at every tested node and ray.
5. **Any grazing issue?**  Batching itself has none.  The audit found and
   corrected a pre-existing in-place mirror/ray-order dependency; forward and
   reversed ray orders are now bitwise identical.
6. **Main spherical bottleneck?**  Legacy `OPMTRX` candidate scanning, not
   `TBINTG_sph`.
7. **Speed-up?**  6.88x to 49.83x in the main real-MARCS cases; 1.35x in the
   deliberately small 983-line NLTE case.
8. **Is NLTE architecture correct?**  Yes.  It uses the unchanged `OPMTRX`
   NLTE opacity/source branch and passed both synthetic and real H-NLTE
   controlled validation.
9. **Safe production combinations?**  Spherical adaptive synthesis with
   `long_continuum=True`, valid precomputed ALMAX/CDR mask and physical ranges,
   and non-legacy mode.
10. **Recommendation?**  Enable in production under the stated gate and keep
    legacy `RKINTS_sph` as compatibility/reference.

## 16. Artifacts

- `analysis/spherical_batched_rkints_benchmark.py`: reproducible benchmark
- `analysis/spherical_batched_rkints_benchmark.json`: scalar results and native timing counters
- `analysis/spherical_batched_rkints_benchmark.npz`: wavelength, per-ray intensity, continuum, and regular-grid flux arrays
