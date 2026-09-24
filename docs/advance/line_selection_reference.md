# Line-Selection Reference

This page is the parameter reference for PySME line selection and line-info
precomputation. It complements [](line_filtering.md), which focuses more on
workflow and examples.

## Mental model

Line selection in PySME is controlled in three layers:

1. `linelist_mode`
   - `all`: synthesize with the full line list
   - `dynamic`: synthesize with a per-segment filtered subset
   - `auto`: deprecated alias of `dynamic`
2. `line_select_method`
   - `internal`: no external CDR/ALMAX metadata path
   - `cdr`: use `central_depth` and `line_range_*`
   - `almax` (default): use `almax_ratio` and `line_range_*`
3. `line_select_policy`
   - `auto`: use method-dependent automatic line-info handling
   - `strict`: require explicit method-specific line-info handling

In practice:

- `linelist_mode` decides whether dynamic filtering is used
- `line_select_method` decides how line metadata is generated/interpreted
- `line_select_policy` decides how strictly that metadata is enforced

## Fixed-grid interval lookup

When valid `line_range_*` metadata is available, SMElib's fixed-grid transfer
path uses a wavelength sweep to visit only the intervals containing the
current wavelength. The candidate indices remain in original line-list order,
so this changes lookup cost rather than opacity summation or numerical line
selection. In automatic mode the sweep is skipped when the ranges are too
broad to save at least roughly 20% of the full scan.

For internally generated fixed-grid metadata, SMElib scans each valid line's
wings until its local line-to-continuum opacity ratio falls below `accrt`.
`GetLineRange` therefore returns physical validity ranges rather than the
temporary `wlcent +/- 150 A` bounds installed when a line list is loaded.

For diagnostic A/B runs, set `SME_INTERVAL_INDEX=0` to disable the sweep or
`SME_INTERVAL_INDEX=1` to force it. Leaving the variable unset (or setting it
to `auto`) uses the automatic cost check.

`ALMAXRange` also leaves its computed line-opacity and Voigt arrays available
for a one-shot hand-off to the next `Transf` call on the same DLL state. The
hand-off requires valid precomputed ranges/masks at the same `accrt`. Updating
the model, abundances, line list, NLTE coefficients, continuum opacity, or
broadening settings invalidates it. This avoids repeating `LINEOPAC` during
first-use precomputation while keeping cached line-info use conservative.

The default non-parallel `linelist_mode="all"` workflow performs a missing or
stale ALMAX calculation inside the main synthesis DLL, after preparing the
first segment's continuum opacity. Its first `Transf` therefore consumes this
one-shot state. Parallel, cache-backed, and dynamic-subsetting workflows use a
separate precompute state and do not receive this first-call reuse.

## Continuum-opacity grid

Continuous opacity is shared by line-info precomputation and final transfer.
By default, PySME evaluates it on an adaptive linear grid with a nominal 1 A
spacing. SMElib inserts known H I thresholds and Mg I/Si I PEACH table knots,
uses an exact 0.02 A guard band around each physical edge, and recursively
splits intervals whose interpolation probes exceed the requested tolerance.
All 13 opacity-source components are cached, so true absorption, coherent
scattering, total extinction, and the scattering source remain consistent.

The user-facing controls are:

- `sme.continuum_grid = "adaptive"` (default): adaptive edge-aware grid
- `sme.continuum_grid = "exact"`: legacy exact `CONTOP` evaluation at every query
- `sme.continuum_grid = 0.5`: fixed edge-aware 0.5 A diagnostic grid
- `sme.continuum_grid_base_step = 1.0`: nominal spacing in A
- `sme.continuum_grid_rtol = 1e-3`: refinement tolerance, measured relative
  to total continuum extinction
- `sme.continuum_grid_min_step = 1e-3`: minimum recursive interval width in A

The refinement test checks true absorption, coherent scattering, and total
extinction, each scaled by total extinction with a small floor. This avoids
refining physically irrelevant components merely because their own value is
close to zero. `"exact"` is intended for reference calculations and numerical
regression tests; normal synthesis does not require choosing a fixed spacing.

## Shared parameters

### `linelist_mode`

Function argument in `solve(...)` and `synthesize_spectrum(...)`.

- `all`: use the full line list
- `dynamic`: use dynamic line filtering
- `auto`: deprecated alias of `dynamic`

### `sme.line_select_method`

- `almax` (default)
- `internal`
- `cdr`

Controls which metadata path is used for line preselection.

### `sme.line_select_policy`

- `auto`
- `strict`

Controls how method-specific line metadata is consumed by the synthesis path.

### `sme.line_select_parallel`

Boolean. Enable or disable parallel metadata updates.

### `sme.line_select_n_jobs`

Worker count for parallel metadata updates.

- `None`: infer automatically
- positive integer: explicit worker count

### `sme.line_select_chunk_size`

Chunk size used when splitting the line list for metadata updates.

This does not affect the default non-parallel, full-line-list ALMAX fast path,
which computes directly in the main synthesis DLL without chunk workers.

### `sme.line_select_recompute`

- `if_stale`: recompute line metadata only when missing or stale
- `always`: recompute every time
- `never`: do not recompute; require existing metadata or cache entries

This is the main control for metadata regeneration.

### `sme.line_select_stale_thres`

Dictionary of stale thresholds, typically including:

- `teff`
- `logg`
- `monh`
- `vmic`
- `accrt`

Used to decide whether previously computed metadata is still valid.

### `sme.line_precompute_database`

Preferred cache-directory parameter.

This is the shared on-disk cache for line precompute products from both `cdr`
and `almax`.

## CDR-specific parameters

### `sme.line_select_cdr_strength_thres`

Threshold for strong-line selection in CDR mode.

### `sme.line_select_cdr_bin_width`

Bin width used in CDR strong-line selection.

CDR currently uses the bin-based strong-line helper directly, so there is no
separate `cdr_use_bins` switch.

## ALMAX-specific parameters

### `sme.line_select_almax_threshold`

Threshold used by ALMAX-based selection.

If `None`, it falls back to `sme.accrt`.

This is a local line-to-continuum opacity-ratio cutoff, not a requested bound
on the final normalized-flux error. Contributions from many individually weak
line wings can accumulate, so the spectrum-level error must be validated for
the intended stellar-parameter and wavelength domain.

### `sme.line_select_almax_use_bins`

Boolean switch controlling which ALMAX strong-line rule is used:

- `False`: simple threshold rule using `almax_ratio >= line_select_almax_threshold`
- `True`: bin-wise cumulative rule using `flag_strong_lines_by_bins(...)`

### `sme.line_select_almax_bin_width`

Bin width used when `line_select_almax_use_bins=True`.

The cumulative bin selector itself lives in SMElib under the CamelCase API
name `SelectStrongLinesByBins`. The Python
`Synthesizer.flag_strong_lines_by_bins` method is retained as a compatibility
wrapper, and both the CDR and binned-ALMAX paths use the native implementation.

## Deprecated or legacy parameters

These are still accepted for backward compatibility, but should not be used in
new code.

### `cdr_database`

Deprecated alias of `line_precompute_database`.

### `sme.cdr_N_line_chunk`

Legacy alias of `sme.line_select_chunk_size`.

### `sme.cdr_parallel`

Legacy alias of `sme.line_select_parallel`.

### `sme.cdr_n_jobs`

Legacy alias of `sme.line_select_n_jobs`.

### `sme.strong_depth_thres`

Legacy alias of `sme.line_select_cdr_strength_thres`.

### `sme.strong_bin_width`

Legacy compatibility field for bin-width based strong-line selection.
Prefer:

- `sme.line_select_cdr_bin_width`
- `sme.line_select_almax_bin_width`

### `sme.line_select_reuse`

Deprecated.

Non-default values currently only trigger a limited internal reuse path by
keeping line opacity around. This is not a fully developed or stable public
cache policy, and new code should leave it at the default `none`.

### `cdr_create`

Legacy-style function argument still used to force regeneration of cached line
metadata products. It remains supported, but is not yet replaced by a clearer
unified name.

## Adaptive transfer-grid semantics

When no `sme.wint` is supplied, plane-parallel synthesis with precomputed
ALMAX or CDR line information constructs the transfer grid in refinement
generations inside native SMElib. Each generation is evaluated through the
indexed fixed-grid opacity path. The initial endpoints, line-centre seeds,
0.3 km/s minimum spacing, and the `accwi` midpoint interpolation criterion
retain the RKINTS definitions.

The active line mask is fixed for the complete transfer calculation. `ALMAX`
or CDR decides whether a line participates, `accrt` defines its wavelength
support, and `accwi` controls wavelength sampling only. In particular,
`accwi` no longer permanently removes a line based on the blended disk-centre
depth at its centre.

Supplying `sme.wint` continues to use that fixed grid directly. Spherical
models and `sme.line_select_method = "internal"` retain legacy RKINTS while
the batched path is validated for those configurations.

`sme.transfer_grid_method = "batched"` is the default. Set it to `"legacy"`
for compatibility or reference calculations. The setting affects only
adaptive transfer; it does not change a supplied fixed `sme.wint` grid.

## Recommended usage

### CDR workflow

```python
sme.line_select_method = "cdr"
sme.line_select_policy = "strict"
sme.line_select_recompute = "if_stale"
sme.line_select_parallel = False
sme.line_select_chunk_size = 2000
sme.line_select_cdr_strength_thres = 0.001
sme.line_select_cdr_bin_width = 0.2
sme.line_precompute_database = "/path/to/cache"
```

### ALMAX workflow

```python
sme.line_select_method = "almax"
sme.line_select_policy = "strict"
sme.line_select_recompute = "if_stale"
sme.line_select_parallel = False
sme.line_select_chunk_size = 2000
sme.line_select_almax_threshold = sme.accrt
sme.line_select_almax_use_bins = False
sme.line_select_almax_bin_width = 0.2
sme.line_precompute_database = "/path/to/cache"
```

## See also

- [](line_filtering.md)
- [](how-to.md)
