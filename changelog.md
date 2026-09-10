# Changelog

## 2026-09-09

- Added `sme.h_stark_convolution = "convolution"` for Br10 and higher.
  The default remains `"legacy"`.
  See the [SMElib v6.13.18 documentation](https://github.com/SpectroscopyMadeEasy/SMElib/blob/v6.13.18/docs/brackett_stark_convolution.md).

## 2026-06-26

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
