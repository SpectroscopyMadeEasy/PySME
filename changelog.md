# Changelog

## v1.2.0 - 2026-09-24

### Performance

- Reworked adaptive synthesis with edge-aware continuum-opacity caching and
  generation-batched, interval-indexed transfer for plane-parallel and
  spherical atmospheres. Complete synthesis was approximately 19--40x faster
  in the canonical 10 A validation matrix; gains vary with line density and
  wavelength coverage.
- Made ALMAX-based line selection the default and reused its immediately
  prepared line-opacity state for the following transfer when safe.
- Reused the resident line list and atmosphere/model during abundance-only
  `solve(...)` iterations while still recomputing abundance-dependent EOS,
  opacity, ALMAX state, and transfer.
- Reduced large-line-list memory through float line-state cache storage with
  double-precision arithmetic and incremental parsing of counted long-format
  VALD `extract stellar` files.

### Correctness

- Removed historical order-dependent second weak-line pruning from the
  optimized path and kept precomputed selected-line masks and physical ranges
  immutable during transfer.
- Corrected fixed-grid physical line ranges, irregular-grid flux integration,
  spherical grazing-ray evaluation order, and contribution-function disk
  integration.
- Invalidated ALMAX metadata on effective abundance changes and recalculated it
  only after the new abundances and ionization state reached SMElib.
- Replaced the historical fixed 400,000-point adaptive-transfer ceiling with
  dynamic capacity sizing for unusually wide, line-rich segments.

### Compatibility

- `sme.transfer_grid_method = "legacy"` retains sequential adaptive transfer
  for compatibility/reference calculations, and
  `sme.line_select_method = "internal"` retains established internal selection.
  User-supplied fixed transfer grids retain the sorted indexed path.
- `sme.continuum_grid = "exact"` retains exact continuum evaluation.
- CDR and optional binned-ALMAX selection now share SMElib's cumulative-bin
  selector, while Python dynamic line-list filtering remains available.
- Clarified that `accrt` is a local line-opacity/support threshold and `accwi`
  is a local adaptive-grid refinement criterion; neither is a global bound on
  final synthesized-spectrum error.

## v1.1.0 - 2026-09-15

### Added

- Added an opt-in continuum-scattering source for plane-parallel and spherical atmospheres; the Planck source remains the default.
- Added the Amarsi & Grevesse (2026) solar abundance pattern as `amarsi2026`.
- Added `sme.nlte.strict` for runs that must stop instead of falling back to LTE when requested NLTE data cannot be applied.
- Added checksum and size validation for downloaded atmosphere and NLTE data.

### Changed

- Interpolate spherical-atmosphere height and radius as a combined logarithmic quantity.
- Build the 3D NLTE hydrogen interpolation grid only when it is first used.
- Pin the bundled SMElib source to v6.13.19.

### Fixed

- Accept VALD-compatible headers that omit the comma after `Vmicro`.
- Fall back to serial CDR line selection when the runtime cannot create worker processes, unless strict line-selection policy is requested.
- Apply progress-bar settings at call time and use the correct Boolean test for disabling progress bars.
- Avoid NumPy shape-assignment warnings in continuum and radial-velocity fitting.
- Close persistence files after failed reads and use a plain Plotly figure outside notebooks.

## v1.0.3 - 2026-09-09

- Added `sme.h_stark_convolution = "convolution"` for Br10 and higher; the default remains `"legacy"`.
- See the [SMElib v6.13.18 documentation](https://github.com/SpectroscopyMadeEasy/SMElib/blob/v6.13.18/docs/brackett_stark_convolution.md).

## 2026-06-26

- Added an opt-in continuum-scattering source treatment for plane-parallel and
  spherical atmospheres using a constant-Eddington-factor moment approximation.
- Fixed the H NLTE abundance-coordinate handling so the standard hydrogen NLTE
  abundance coordinate remains stable during synthesis.
- Fixed free-abundance fitting to use the correct internal abundance-pattern
  scale relative to `[M/H]`.
- Fixed spherical atmosphere interpolation so `height` is interpolated
  consistently with the other atmospheric structure quantities.
- Fixed derived abundance parameter handling in `solve()` so abundance keys are
  parsed consistently, including capitalized forms such as `"Abund Ti"`.
- Improved SMElib robustness:
  - tolerate missing HLINOP warning symbols in older SMElib builds
  - rebuild the DLL object inside multiprocessing worker processes when needed
- Added explicit abundance-scale views and updated documentation around
  `sme.abund.A[...]` and `sme.abund.pattern[...]`.
- Unified line-selection controls around the `line_select_*` interface while
  keeping legacy `cdr_*` pathways available for compatibility.
- Added an experimental `profile_nlte` interface for profile-based NLTE
  corrections. This path is disabled by default and remains experimental.
- Compatibility and deprecations:
  - `dynamic_param` is deprecated in favor of `derived_param`
  - `cdr_database` is deprecated in favor of `line_precompute_database`
  - `linelist_mode="auto"` is deprecated in favor of `linelist_mode="dynamic"`
  - direct abundance assignment through `sme.abund["X"]` is still supported,
    but now emits a warning because it writes the internal pattern rather than
    the final abundance used in synthesis
 
## 2026-02-10

- Added `derived_param` as the preferred name for derived-parameter callbacks in `solve()`.
- Kept backward compatibility with `dynamic_param`:
  - If `dynamic_param` is used, a `DeprecationWarning` is emitted.
  - If both `derived_param` and `dynamic_param` are provided (and differ), `ValueError` is raised.
- Updated internal solver logic and user-facing messages to use the "derived parameter" terminology.
- Added `smelib_lineinfo_mode` passthrough in `solve()` call paths (`_residuals` and `_jacobian`) so fitting runs can use SMElib precomputed line-info modes.
- Updated `linelist_mode` naming:
  - Preferred values are now `"all"` and `"dynamic"`.
  - `"auto"` is kept as a deprecated compatibility alias and maps to `"dynamic"` with a `DeprecationWarning`.
- Added segment-aware optional input `sme.wint` for synthesis transfer grids.
  - Priority is now `sme.wint[segment]` first, then internal cached grids (when enabled), then SMElib adaptive grid generation.
- Updated user docs accordingly (`sme_struct`, `quickstart`, and `how-to`).
