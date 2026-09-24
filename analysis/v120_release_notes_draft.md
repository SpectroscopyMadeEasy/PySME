# PySME v1.2.0 - faster synthesis and improved numerical robustness

> Draft release notes. Do not publish from this file.

## Highlights

- Complete synthesis was approximately 19--40x faster in the canonical 10 A
  validation matrix; gains depend on line density, wavelength coverage, and
  the previous bottleneck.
- The optimized adaptive transfer path supports plane-parallel and spherical
  atmospheres and retains existing LTE, tested NLTE, and hydrogen-profile
  physics.
- Continuum-opacity caching and interval-indexed transfer remove repeated work
  without replacing PySME's adaptive internal wavelength grid.
- Float line-state caches and incremental parsing of supported large VALD files
  reduce memory pressure for very large line lists.
- v1.2 corrects historical weak-line, line-range, irregular-grid integration,
  and spherical ray-order behavior.

## Upgrading from v1.1.x

Standard synthesis scripts generally require no changes. The default ALMAX
workflow computes missing or stale line information and automatically uses the
optimized path when its mask and physical ranges are valid.

`sme.wave` is the requested output/observation grid and does not force fixed
native transfer. Set `sme.wint` only when the transfer calculation itself must
use a supplied fixed grid.

| Existing workflow | Code changes required? | v1.2 behavior | Can results differ from v1.1.x? |
|---|---|---|---|
| Default adaptive ALMAX synthesis | No | Optimized automatically | Possibly slightly, due to numerical and correctness changes |
| Adaptive synthesis with valid CDR metadata | No | Optimized automatically | Possibly, especially near a cumulative selection boundary |
| Only `sme.wave` supplied | No | Still optimized adaptive transfer | Same as the selected ALMAX/CDR workflow |
| Explicit `sme.wint` | No | Sorted, indexed fixed-grid transfer | Only from applicable v1.2 changes |
| `line_select_method="internal"` | No | Established internal/legacy adaptive path | Other v1.2 fixes still apply |
| `transfer_grid_method="legacy"` | No | Sequential legacy adaptive scheduler | Compatibility mode; not a complete v1.1 rollback |

For historical comparisons, the supported controls are:

```python
sme.transfer_grid_method = "legacy"
sme.line_select_method = "internal"
sme.continuum_grid = "exact"
```

Use only the controls relevant to the comparison. They do not undo the
irregular-grid integration fix, spherical ray-order fix, or unrelated bug
fixes, so bitwise v1.1.x reproduction is not guaranteed.

## Performance

The canonical complete-synthesis validation improved by approximately
19--40x over the historical path for 10 A cases. Tested Mg b and H-alpha cases
improved by approximately 18--43x. Sparse or compact workloads can gain less;
wide line-rich spectra can see larger absolute time savings.

### Adaptive continuum-opacity cache

The default edge-aware adaptive cache replaces repeated exact continuum
queries with interpolation over refined opacity nodes. In representative
release cases, hundreds of continuum queries were served per exact evaluation.
The standalone solar `ALMAXRange + CentralDepth` continuum/CDR task improved
from 61.23 s to 3.56 s, or about 17.2x.

In the controlled direct-synthesis matrix, the maximum normalized-flux change
from continuum interpolation was about `1.4e-9`.

### Generation-batched adaptive transfer

The scheduler evaluates sorted refinement generations through physical
line-range indexing. It does not remove the physical transfer-wavelength
evaluations: for example, the cool-dwarf comparison evaluated the same 10,049
wavelengths in both schedulers. Candidate-line visits decreased by
approximately 99.63--99.94% across the release matrix.

Representative native transfer speed-ups were 23.9x for the dense
plane-parallel cool-dwarf case and 6.9--31x across the principal spherical
cases, with larger transfer-only gains in some Mg b and H-alpha tests.
Controlled outputs were bitwise identical.

The continuum-cache and batching component factors are not multiplicative;
they remove different bottlenecks that dominate different stellar models.

## Correctness changes

- The optimized path no longer applies the historical sequential `MARK=2`
  weak-line pruning after ALMAX/CDR selection. In a tested metal-poor case, the
  old second pruning moved the spectrum away from a minimally pruned reference.
- Precomputed selected-line masks and physical wavelength ranges remain
  immutable during optimized transfer.
- Angle-dependent line and continuum intensities are resampled to a common
  regular log-wavelength grid before disk integration and velocity broadening.
- Spherical normal and grazing rays are evaluated in a deterministic safe order,
  removing dependence on caller `mu` ordering.
- Internally generated transfer capacity is sized dynamically when the
  historical 400,000-point allocation is insufficient.

Small differences from older releases in affected spectra are intentional
correctness changes, not accuracy traded for speed. Not every spectrum is
affected.

## Numerical validation

- Controlled batching comparisons were bitwise identical for complete spectra
  and native common-node `sint`/`cint`, including plane-parallel, spherical
  grazing-ray, hydrogen-profile, and tested H-NLTE cases.
- Adaptive continuum interpolation changed direct-synthesis normalized flux by
  at most about `1.4e-9` in the release matrix.
- A separate cumulative-bin CDR selection audit reached about `7.24e-6` at a
  weak-line ranking boundary; this was a selection decision, not transfer
  error.
- Float line-state storage, with double-precision arithmetic, changed tested
  normalized flux by at most about `1.6e-8` and changed no selected-line masks
  or physical ranges.
- The release source suite passed 188 tests with one expected failure in the
  recorded validation environment.

## Memory

`LINEOP`, `AVOIGT`, and `VVOIGT` cache storage changed from double to float,
halving their logical payload. In the 1.26-million-line, 55-depth validation
case this reduced the three-cache payload from about 1.66 GB to 0.83 GB and
reduced synthesis-stage RSS by about 0.81 GB.

Counted VALD `long + extract stellar` files are now parsed incrementally. The
same 1.26-million-line case reduced parser peak RSS from about 4.19 GB to
1.86 GB, or about 56%, while preserving the returned data frame.

## Compatibility and fallbacks

- `sme.transfer_grid_method="legacy"` retains sequential adaptive transfer for
  compatibility and reference calculations.
- `line_select_method="internal"` retains the established internal selection
  and adaptive path.
- A supplied `sme.wint` uses fixed-grid indexed transfer; adaptive batching is
  unnecessary because the complete sorted wavelength grid is known.
- High-level synthesis normally recomputes missing or stale ALMAX/CDR metadata.
  Automatic policy can fall back if valid metadata cannot be supplied; strict
  policy reports the problem.
- The low-level `long_continuum=False` combination retains the established
  path and is not a high-level `SME_Structure` configuration.

## Known limitations

- SMElib remains process-global and non-reentrant. High-level threaded calls
  are serialized for correctness; use processes for true parallel synthesis.
- Very wide, million-line, single-segment synthesis can still require several
  GiB per process.
- Generation batching is not applied to internal line selection or explicit
  legacy adaptive transfer.
- Legacy controls preserve selected historical semantics but do not undo every
  v1.2 correctness fix or guarantee bitwise v1.1.x output.
