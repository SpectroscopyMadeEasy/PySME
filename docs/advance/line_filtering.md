# Line Filtering

For wide wavelength coverage (or many segments), using the full line list in
every segment is expensive. PySME provides dynamic line filtering to keep only
relevant lines per segment (see [Jian et al. (2026)](https://ui.adsabs.harvard.edu/abs/2026arXiv260504007J/abstract)).

## Core Options

Use `linelist_mode` and `line_select_method` together in synthesis or solve:

- `"all"`: use all lines (default).
- `"dynamic"`: filter lines by precomputed line properties (recommended for long spectra).
- `"auto"`: legacy alias of `"dynamic"` (deprecated).
- `line_select_method="almax"`: use `almax_ratio` + `line_range_*` (default).
- `line_select_method="internal"`: use SMElib's legacy internal selection.
- `line_select_method="cdr"`: use `central_depth` + `line_range_*`.

## How Dynamic Filtering Works

When `linelist_mode="dynamic"`:

1. PySME ensures line metadata exists for the selected method.
   - `cdr`: `central_depth`, `line_range_s`, `line_range_e`
   - `almax`: `almax_ratio`, `line_range_s`, `line_range_e`
2. If metadata is missing or stale, PySME updates it automatically.
2. For each segment, PySME keeps lines that overlap the segment range
   (with broadening margin) and pass a strength threshold.
3. Only this reduced line subset is sent to SMElib for that segment.

This can significantly reduce runtime for long or segmented spectra.

For a complete parameter-by-parameter reference, including deprecated aliases
and recommended replacements, see [](line_selection_reference.md).

## Main Controls

- Shared controls:
  - `sme.line_select_method`: `internal | cdr | almax`
  - `sme.line_select_policy`: `auto | strict`
  - `sme.line_select_parallel`: enable/disable parallel metadata update
  - `sme.line_select_n_jobs`: number of worker processes
  - `sme.line_select_chunk_size`: chunk size for metadata update
  - `sme.line_select_recompute`: `if_stale | always | never`
  - `sme.line_select_stale_thres`: stale thresholds for model changes
- CDR method controls:
  - `sme.line_select_cdr_strength_thres`
  - `sme.line_select_cdr_bin_width`
  - `line_precompute_database` / `cdr_create` (function args)
- ALMAX method controls:
  - `sme.line_select_almax_threshold`
  - `sme.line_select_almax_use_bins`
  - `sme.line_select_almax_bin_width`

`line_precompute_database` is a shared on-disk cache for both `cdr` and
`almax`. Cache entries are isolated by `(method, linelist_hash, stellar params)`,
so one folder can safely store multiple linelists and both methods together.
Legacy `cdr_database` is still accepted as a deprecated alias.

The current on-disk cache key does not encode a custom element-by-element
abundance pattern. After ALMAX metadata has been associated with an in-memory
line list, PySME detects an abundance change and bypasses that cache. A cache
built under a different custom abundance pattern should not be supplied on the
first synthesis; regenerate it or disable `line_precompute_database`.

### Recompute vs. reuse

- `line_select_recompute` controls whether line metadata is recomputed when it
  is missing or stale:
  - `if_stale`: recompute only when needed
  - `always`: always recompute
  - `never`: require existing metadata or cache entries
- ALMAX staleness includes an exact comparison of the effective elemental
  abundance vector. Changing any abundance therefore triggers new ALMAX ratios,
  strong-line flags, and validity ranges. This abundance check is independent
  of the approximate atmosphere thresholds in `line_select_stale_thres`.
- `line_select_reuse` is deprecated. Non-default values still enable a limited
  internal reuse path by keeping line opacity around, but this is not a fully
  developed cache policy and should not be treated as a stable public API.

### ALMAX strong-line rule

- `line_select_almax_use_bins=False`:
  - `strong = (almax_ratio >= line_select_almax_threshold)`
- `line_select_almax_use_bins=True`:
  - `strong = flag_strong_lines_by_bins(wl, almax_ratio, threshold=line_select_almax_threshold, bin_width=line_select_almax_bin_width)`

`line_select_almax_threshold` is the single ALMAX threshold parameter for both
rules. If it is `None`, it falls back to `sme.accrt` (legacy-compatible
behavior).

`accrt` and `line_select_almax_threshold` are local line-to-continuum opacity
ratio cutoffs, not bounds on the final normalized-flux error. The default is
`1e-4`; accumulated weak-line contributions can produce a larger flux change.

For a missing or stale ALMAX result in non-parallel `"all"` mode, PySME runs
`ALMAXRange` in the synthesis DLL and immediately reuses its line-opacity and
Voigt state in the first transfer calculation. Cached, parallel, and
`"dynamic"` workflows retain their separate precompute path.

During an abundance-only solve, PySME may retain the resident line list and
atmosphere, but it does not retain abundance-dependent ALMAX results. The new
abundances and ionization state are installed first, then `ALMAXRange` is
rerun. Prepared-state reuse is not enabled for CDR selection because CDR cache
metadata is not currently abundance-aware.

CDR does not have a separate `use_bins` switch because its current strong-line
selection already uses the bin-based helper internally.

The cumulative bin rule is implemented once in SMElib as
`SelectStrongLinesByBins`. Both CDR and optional binned-ALMAX selection call
that native implementation. With `linelist_mode="dynamic"`, Python still uses
the returned mask and ranges to reduce the line list before passing it to the
main synthesis DLL; parallel CDR calculation and cache handling are unchanged.
This first-stage integration does not change the default individual ALMAX rule
or solve its accumulated-weak-line limitation.

## Example 1: Dynamic Filtering in Synthesis (CDR)

```py
from pysme.synthesize import Synthesizer, synthesize_spectrum

synth = Synthesizer()
sme = synth.update_cdr(sme)              # populate central_depth / line_range_* once
sme.line_select_method = "cdr"
sme.line_select_cdr_strength_thres = 0.02

sme = synthesize_spectrum(sme, linelist_mode="dynamic")
```

## Example 2: Dynamic Filtering in Synthesis (ALMAX)

```py
from pysme.synthesize import synthesize_spectrum

sme.line_select_method = "almax"
sme.line_select_almax_threshold = sme.accrt
sme.line_select_almax_use_bins = True
sme.line_select_almax_bin_width = 0.2

sme = synthesize_spectrum(sme, linelist_mode="dynamic")
```

## Example 3: Dynamic Filtering in Solve

```py
from pysme.solve import solve

fit = ["teff", "logg", "monh", "vmic"]
sme.line_select_method = "almax"
sme.line_select_almax_use_bins = True

sme = solve(
    sme,
    fit,
    linelist_mode="dynamic",
    line_precompute_database="path/to/line_precompute_db",
    cdr_create=False,                     # set True to force regeneration
)
```

## Practical Guidance

- The default is `line_select_method="almax"` with
  `line_select_almax_threshold = None`, which uses `sme.accrt`.
- Set `line_select_method="internal"` to reproduce the legacy selection path.
- Enable `line_select_almax_use_bins=True` when you want bin-wise cumulative pruning.
- Use `"all"` for short, narrow windows where filtering overhead may not help.
- Use `"dynamic"` for wide ranges or many segments.
