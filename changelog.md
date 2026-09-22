# Changelog

## Unreleased

### Changed

- Reused the resident SMElib line list and atmosphere/model during
  abundance-only `solve(...)` iterations. Abundances, EOS, opacity, and
  transfer are still recomputed for every trial; mixed-parameter fits and
  workflows requiring line-selection recomputation retain the full setup path.
- Made the adaptive continuum-opacity grid the default for line-info
  precomputation and synthesis. It uses a nominal 1 A grid, physical
  H I/Mg I/Si I knots, exact edge guard bands, and recursive curvature
  refinement. `sme.continuum_grid = "exact"` retains the reference path.
- Fixed internal fixed-grid transfer so `GetLineRange` returns opacity-based
  validity ranges instead of the `wlcent +/- 150 A` placeholders initialized
  by `InputLineList`. This makes CDR range metadata usable by interval-based
  line selection.
- Made ALMAX-based line selection the default. A missing or stale ALMAX result
  in the non-parallel, full-line-list workflow is now calculated in the main
  SMElib instance so the first transfer reuses its line-opacity and Voigt
  state. Set `sme.line_select_method = "internal"` for the legacy behavior.
- Clarified that `accrt` is a local line-to-continuum opacity-ratio threshold,
  not a bound on the final synthesized-spectrum error.
- Moved cumulative wavelength-bin line selection into SMElib as
  `SelectStrongLinesByBins`. CDR and optional binned-ALMAX selection now share
  this native implementation while Python dynamic line-list pruning remains
  available.

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
