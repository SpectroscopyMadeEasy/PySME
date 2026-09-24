# Generation-batched RKINTS prototype

## Executive conclusion

The performance hypothesis is confirmed, with one important qualification.

For the 5195--5205 A cool metal-rich case, the legacy adaptive transfer
calculation took **118.70 s**, while the generation-batched prototype took
**4.91 s** for the same 10,049 evaluated wavelengths: a **24.17x transfer
speed-up**.  The final wavelength arrays and normalized intrinsic fluxes were
identical at every node.

Native profiling shows that this is not mainly a reduction in radiative-
transfer evaluations.  Both paths made approximately 10,050 `OPMTRX` calls.
The difference is that the pre-set-grid path can use `LineIntervalSweep` and
`OPMTRXIndexed`: it reduced the candidate line visits by 99.66%, cutting
`OPMTRX` from 114.19 s to 3.26 s.  `TBINTG` was already cheap and changed only
from 0.064 s to 0.055 s.

The follow-up second-pruning audit resolved the semantic qualification:
ALMAX-only transfer is substantially closer to the minimally pruned reference,
whereas legacy `MARK=2` removal is order-dependent. The generation-batched
scheduler now runs directly inside production SMElib for plane-parallel
adaptive synthesis with precomputed ALMAX/CDR line information. Its line mask
and physical ranges remain immutable.

The production implementation reproduces the prototype wavelengths and fluxes
exactly in all four models. Fixed grids, spherical models, and legacy internal
line selection retain their established paths.

## Scope and benchmark setup

- Branch: `codex/cdr-physical-range-selection`
- Window returned to the user: 5195--5205 A
- SMElib opacity window: approximately 5192.48--5207.52 A
- Line list: `/tmp/pysme_line_centric_windows/window_0010A_pad_150A.lin`
  (378,154 input lines)
- `accrt = 1e-4`
- `accwi = 3e-3`
- RKINTS minimum velocity spacing: 0.3 km/s
- Seven limb angles from the normal PySME setup
- `vmac = vsini = 0`; instrumental resolution disabled for the intrinsic
  transfer calculation
- Accuracy reference: a 0.000625 A fixed grid (16,001 points in the returned
  10 A window)
- Timings are single-process wall times on the same machine and library state.

The original experiment remains reproducible in
`analysis/batched_rkints_benchmark.py`. The production scheduler now lives in
native `Transf`/`RKINTS_batched`: each generation uses the existing indexed
opacity and transfer routines without duplicating opacity or RT physics.

## Legacy execution audit

The plane-parallel adaptive path currently executes in this order:

1. `Transf` restores the precomputed ALMAX strong mask and physical line
   ranges, then computes `LINEOPAC` for retained lines.
2. RKINTS evaluates the blue endpoint.
3. It scans the SMElib line list in its existing order.  For every active line
   centre separated by more than 0.3 km/s from the last grid point, it inserts
   the preceding midpoint and the line centre and immediately evaluates both.
4. At the line centre it applies
   `1 - I_line(mu=1) / I_cont(mu=1) < accwi`; if true, it sets that line's
   `MARK=2`.  Later redder evaluations therefore see a changed active-line
   state.
5. It adds the red endpoint.
6. It walks intervals left-to-right.  Every interval midpoint is evaluated and
   retained.  The interval is accepted when

   ```text
   (|I_mid - (I_left + I_right)/2|
    + 0.005 |I_left - I_right|) / I_cont_mid < accwi
   ```

   or when its left half-width is no larger than the 0.3 km/s minimum.
   Rejected intervals are refined depth-first.

The disk-centre ray is selected by the largest `mu`, independent of caller
ordering.  With `long_continuum=True`, the midpoint's exact continuum is used
as the normalization.

The pre-set-grid branch instead sees the complete sorted wavelength array in
advance.  It builds a `LineIntervalSweep` from the physical validity ranges,
then calls `OPMTRXIndexed` at each wavelength.  Opacity is still accumulated in
original line-list order; only lines whose ranges overlap the current
wavelength are visited.

## Prototype algorithm

The diagnostic scheduler retains the legacy seed and refinement geometry:

1. Generate the blue endpoint plus the same `(previous + line centre)/2`, line
   centre pairs, followed by the red endpoint.
2. Evaluate that complete initial grid with one fixed-grid `Transf` call.
3. Store every returned line and continuum intensity by exact wavelength.
4. For every active interval, generate its midpoint.  Evaluate all new sorted
   midpoints in one fixed-grid `Transf` call with `keep_lineop=True`.
5. Apply the exact legacy disk-centre error expression and minimum spacing.
6. Retain every tested midpoint.  Rejected intervals create the next
   generation's left and right children.
7. Reuse cached endpoints; only newly generated midpoint wavelengths are sent
   to the next `Transf` call.

The sole intentional semantic difference is that the strict precomputed ALMAX
mask remains fixed.  The prototype does not reproduce the sequential
line-centre `MARK=2` mutation.

## Runtime results

| model | legacy points | batch points | midpoint generations | fixed-grid batches | legacy transfer | batch transfer | transfer speed-up |
|---|---:|---:|---:|---:|---:|---:|---:|
| Solar dwarf (5777/4.44/0.0) | 2,447 | 2,447 | 4 | 5 | 1.503 s | 0.147 s | 10.20x |
| Metal-poor dwarf (6000/4.0/-2.0) | 1,245 | 1,243 | 6 | 7 | 0.331 s | 0.057 s | 5.78x |
| Cool metal-rich dwarf (4250/4.5/+0.3) | 10,049 | 10,049 | 1 | 2 | 118.704 s | 4.912 s | 24.17x |
| Giant (4500/2.0/0.0) | 9,933 | 9,933 | 1 | 2 | 102.437 s | 4.375 s | 23.41x |

The shared PySME model/ALMAX setup took 1.96--3.08 s.  Adding this measured
setup to each transfer time gives an indicative (not production end-to-end)
speed-up of 1.64x solar, 1.14x metal-poor, 15.25x cool metal-rich, and 14.63x
giant.  This sum includes the benchmark's 201-point initialization transfer,
so it should not be presented as a finalized production `synthesize_spectrum`
timing.

Python grid bookkeeping was insignificant: 0.0016 s solar, 0.0008 s
metal-poor, 0.0050 s cool dwarf, and 0.0057 s giant.

## Refinement generations

The evaluated wavelength counts, including the initial grid as generation 0,
were:

| model | evaluated wavelengths per generation |
|---|---|
| Solar dwarf | 952, 951, 252, 190, 102 |
| Metal-poor dwarf | 488, 487, 80, 96, 66, 18, 8 |
| Cool metal-rich dwarf | 5,025, 5,024 |
| Giant | 4,967, 4,966 |

The apparently non-monotonic 80 -> 96 step is expected: every rejected parent
creates two children, while many other parents are accepted in the same
generation.

No cached endpoint was re-evaluated.  The number of unique wavelength
evaluations equals the final batch grid size in every model.

## Native cool-star profile

`SME_TIMING=1` gave the following decomposition for the pathological cool
case:

| component | legacy adaptive | batched initial + midpoint | observation |
|---|---:|---:|---|
| Total `Transf` | 115.734 s | 4.836 s | 23.9x |
| Setup / `LINEOPAC` | 1.414 / 1.401 s | 1.429 / 1.416 s | unchanged |
| RKINTS body | 114.319 s | 1.708 + 1.698 s | interval-index path used by batches |
| `OPMTRX` | 114.194 s, 10,050 calls | 1.635 + 1.627 s, 10,051 calls | same call count, far fewer candidate lines |
| `TBINTG` | 0.064 s, 20,100 calls | 0.028 + 0.028 s, 20,100 calls | not the bottleneck |
| Interval candidates | full scan | 8.59 million vs 2.524 billion full-scan visits | 99.66% reduction |

Thus “batching” is valuable because it exposes a sorted generation grid to the
existing interval index.  Merely reducing C/Python call overhead would not
explain the result.

## Grid equivalence

| model | common nodes | legacy-only | batch-only | Hausdorff distance |
|---|---:|---:|---:|---:|
| Solar dwarf | 2,447 | 0 | 0 | 0 A |
| Metal-poor dwarf | 1,243 | 2 | 0 | 0.006794 A |
| Cool metal-rich dwarf | 10,049 | 0 | 0 | 0 A |
| Giant | 9,933 | 0 | 0 | 0 A |

The two extra metal-poor legacy nodes arise from a refinement decision changing
after the sequential weak-line state diverges.  They are not a failure to
reproduce the midpoint construction itself.

## Batch versus legacy spectrum

Both spectra were interpolated to the same 0.000625 A regular grid before
comparison and convolution.  Intrinsic flux is a direct projected-area sum of
the limb intensities; no irregular grid was passed to `integrate_flux`.

| model | max intrinsic delta | intrinsic RMS | max delta R=20k | max delta R=60k | max absolute EW delta |
|---|---:|---:|---:|---:|---:|
| Solar dwarf | 1.85e-5 | 1.31e-6 | 4.51e-6 | 1.02e-5 | 2.00e-8 A |
| Metal-poor dwarf | 1.46e-3 | 1.07e-4 | 2.79e-4 | 5.99e-4 | 2.30e-5 A |
| Cool metal-rich dwarf | 0 | 0 | 0 | 0 | 0 A |
| Giant | 1.98e-7 | 5.33e-9 | 8.94e-9 | 2.56e-8 | 1.20e-9 A |

At common native nodes the maximum batch-minus-legacy differences were zero
for the cool dwarf, `1.98e-7` for the giant, `1.86e-5` for the Sun, and
`2.71e-3` for the metal-poor dwarf.  The larger common-node result in the last
case proves that the discrepancy is line-state history, not interpolation or
the two missing refinement nodes.

## Adaptive grid versus fixed-grid reference

This comparison measures the existing `accwi=3e-3` adaptive-grid interpolation
accuracy, not a new batch approximation.  Where the batch and legacy grids are
identical, their values are identical or nearly so.

| model | max intrinsic error | intrinsic RMS | max error R=20k | max error R=60k |
|---|---:|---:|---:|---:|
| Solar dwarf | 5.98e-3 | 4.40e-4 | 1.04e-3 | 9.95e-4 |
| Metal-poor dwarf | 4.91e-3 | 3.18e-4 | 4.69e-4 | 1.09e-3 |
| Cool metal-rich dwarf | 8.36e-4 | 5.92e-5 | 2.69e-5 | 6.28e-5 |
| Giant | 1.03e-3 | 6.43e-5 | 2.68e-5 | 6.65e-5 |

This reiterates that `accwi` is a local disk-centre midpoint heuristic, not a
global bound on the final interpolated spectrum.

## Production decision

The prototype met the performance target and proved that the existing
fixed-grid interval index can remove the cool/giant pathology without reducing
the number of physical transfer evaluations.  It also reproduces the exact
legacy grid in three of four models.

The follow-up accuracy audit selected deterministic ALMAX semantics: the
second pruning moved the metal-poor result away from a minimally pruned
reference and could change with decision order. Production adaptive transfer
therefore keeps the precomputed mask and ranges fixed and uses `accwi` only for
wavelength refinement.

The final native production validation reproduced all prototype wavelengths
and normalized fluxes exactly, and preserved every physical range bitwise.
Measured native/legacy transfer times were 0.148/1.512 s solar, 0.055/0.320 s
metal-poor, 4.818/115.229 s cool metal-rich, and 4.351/102.385 s giant. These
correspond to 10.24x, 5.85x, 23.91x, and 23.53x transfer speed-ups. The legacy
path remains available through `sme.transfer_grid_method = "legacy"` and is
also retained automatically for internal line selection and spherical
transfer.

## Artifacts

- `analysis/batched_rkints_benchmark.py`: reproducible prototype and benchmark
- `analysis/batched_rkints_benchmark.json`: scalar results and generation logs
- `analysis/batched_rkints_benchmark.npz`: wavelength/flux arrays
- `analysis/batched_rkints_cool_timing.log`: native cool-star timing breakdown
- `analysis/production_batched_rkints_validation.json`: production/prototype
  exact-equality check
