# PySME v1.2.0 documentation audit

## 1. Scope and files inspected

This audit documents the v1.2 release-candidate behavior without changing
scientific code. The source-of-truth inspection covered:

- `src/pysme/sme.py`: public structure fields and defaults;
- `src/pysme/synthesize.py`: line-selection configuration, metadata refresh,
  production path setup, `wave`/`wint` priority, intensity handling, and flux
  integration;
- `src/pysme/solve.py`: solve defaults, state reuse, and high-level locking;
- `src/pysme/sme_synth.py` and `src/pysme/sme_synth_cw.py`: current and legacy
  low-level wrapper behavior;
- `smelib/src/sme/sme_synth_faster.cpp`: batched gating, immutable line state,
  continuum cache, geometry-specific transfer, dynamic storage, and float
  cache declarations;
- `test/test_dll.py`, `test/test_smelib_thread_safety.py`,
  `test/test_halpha_regression.py`, and `test/test_sme_structure.py`;
- all existing user/developer pages containing the audited controls;
- `analysis/v120_release_performance_summary.md`,
  `analysis/batched_rkints_benchmark.md`,
  `analysis/spherical_batched_rkints_benchmark.md`,
  `analysis/v120_line_state_memory_audit.md`, and
  `analysis/v120_vald_parser_memory_audit.md`.

## 2. Source-of-truth API and default audit

| Concept | Current source/API | Current docs location | Outdated/missing before audit? | Action |
|---|---|---|---|---|
| `line_select_method` | `SME_Structure` attribute; `almax` default; accepts `almax`, `cdr`, `internal` | Line Filtering; Line-Selection Reference | Present, but not tied clearly to support/sampling | Added canonical three-stage model and optimized-path table |
| `linelist_mode` | Function argument to synthesis/solve; `all` default; `dynamic`; deprecated `auto` alias | Line Filtering; Line-Selection Reference | Present | Clarified that both `all` and `dynamic` can use batched transfer |
| `line_select_policy` | `SME_Structure` attribute; `auto` default or `strict` | Line-Selection Reference | Present but fallback consequence was not summarized | Added automatic recompute/fallback versus strict-error behavior |
| `accrt` | Public field; default `1e-4` | SME Structure; line-selection pages | Mostly current | Standardized as a local line-opacity/support threshold, not a final-flux bound |
| `accwi` | Public field; default `3e-3` | SME Structure; Line-Selection Reference | Mostly current | Standardized as a local emergent-intensity refinement criterion; no optimized second pruning |
| `long_continuum` | Low-level `SME_DLL.Transf` keyword; default `True`; not an `SME_Structure` field | Low-level autodoc only | Could have been mistaken for a user configuration | Documented only as a developer fallback condition |
| `transfer_grid_method` | Public field; default `batched`; accepts `batched`, `legacy` | SME Structure; Line-Selection Reference | SME Structure incorrectly said plane-parallel only | Corrected to plane-parallel and spherical; documented legacy compatibility behavior |
| Legacy adaptive transfer | Publicly selected with `sme.transfer_grid_method="legacy"` | Line-Selection Reference | Present but migration meaning missing | Added supported reproduction guidance and limitations |
| ALMAX/CDR mask/ranges | Python columns `strong`, `line_range_s/e`; passed to native precomputed line info; normally regenerated if stale | Line Filtering; Line-Selection Reference | Implementation/immutability split was scattered | Documented immutable optimized semantics and refresh behavior |
| `specific_intensities_only` | Public field; default `False` | Flux and intensity | Present | Added native-grid/resampling and skipped post-processing details |
| `wint` | Public optional fixed native transfer grid; `None` by default; also intensity output when requested | SME Structure; Flux and intensity | Incorrectly described as an “adaptive grid” input | Corrected fixed-grid meaning and distinguished it from `wave` |
| `sint`, `cint` | Public segmented native intensity outputs, normally populated for `specific_intensities_only=True` | Flux and intensity | Present but internal-grid semantics incomplete | Added shapes, grid relationship, and post-processing behavior |
| `wave` versus `wint` | `wave` is requested output/observation grid; only `wint` is passed as a fixed native grid | Previously implicit | Missing and user-confusing | Added to migration and permanent path tables |
| Continuum-grid controls | Public attributes: `adaptive`, `exact`, or positive fixed step; defaults 1 A, `rtol=1e-3`, minimum 0.001 A | Line-Selection Reference | Present only on advanced line-selection page | Added permanent performance page and SME Structure entries |
| `SME_DLL.session()` | Exists on the low-level wrapper as a transaction lock; high-level workflows use a process-wide decorator | Not documented as user API | Public intent is not established | Deliberately not promoted; documented only the supported high-level concurrency behavior |

The exact production batched gate is: internally generated transfer grid,
long continuum, `transfer_grid_method="batched"`, and validated precomputed
mask/ranges. It applies to both plane-parallel and spherical atmospheres.

## 3. Outdated or missing documentation found

- There was no v1.2 migration page answering whether existing code changes,
  which paths accelerate automatically, why results can change, or what legacy
  controls do.
- Existing docs did not state clearly that supplying only `sme.wave` still uses
  adaptive native transfer.
- `sme.wint` was called an “adaptive grid” input even though supplying it fixes
  the native transfer wavelengths.
- `transfer_grid_method` was described as plane-parallel only after spherical
  batching had become production behavior.
- The relationship between line membership, physical support, and wavelength
  sampling was spread over several pages.
- The flux page did not describe common regular log-wavelength resampling before
  disk integration/broadening.
- The old developer note said transfer range discovery scans all input lines
  without describing the precomputed/indexed production path.
- There was no permanent optimized-path overview, synthesis architecture page,
  large-VALD memory note, or parallelism FAQ.
- The changelog did not yet include float line-state storage and was organized
  as a long implementation-detail list rather than performance/correctness/
  compatibility guidance.

## 4. Files changed

User-facing documentation:

- `docs/index.rst`
- `docs/getting_started/index.rst`
- `docs/getting_started/whats_new_v120.md`
- `docs/advance/index.rst`
- `docs/advance/synthesis_performance.md`
- `docs/advance/line_filtering.md`
- `docs/advance/line_selection_reference.md`
- `docs/advance/faq.md`
- `docs/advance/fordev.md`
- `docs/fundamentals/flux_inten.md`
- `docs/fundamentals/linelist.md`
- `docs/concepts/sme_struct.md`
- `README.md`

Developer/release documentation:

- `docs/dev/index.rst`
- `docs/dev/synthesis_engine.md`
- `docs/dev/changelog.md`
- `changelog.md`
- `analysis/v120_release_notes_draft.md`
- this audit report.

## 5. Release overview and migration guidance added

The new “PySME v1.2” page is the release landing page for all readers. It
summarizes user-visible value, action required, scientific-result expectations,
and the two conceptual sources of the speed-up without introducing synthesis
implementation vocabulary. It then routes existing users, advanced users, and
developers to progressively more detailed pages.

Essential migration guidance is included directly in the overview: standard
scripts generally need no changes, LTE analyses do not generally need rerunning
solely because of the performance work, and the separate historical NLTE
notice still applies. The path-selection matrix, exact controls, detailed
validation, and implementation methodology remain on the advanced and
developer pages instead of being duplicated in a release-specific migration
page.

The homepage links directly to this page and retains the existing NLTE
correctness warning.

## 6. Permanent behavior documentation added

The new “Synthesis performance and numerical controls” page contains the
current path-selection table, user-level pipeline, continuum-cache behavior,
generation batching, sparse line evaluation, canonical `ALMAX/accrt/accwi`
semantics, and conservative performance expectations. It is version-agnostic
current behavior rather than a migration ledger.

Line-selection, flux/intensity, SME Structure, VALD, and FAQ pages were updated
at their existing points of authority rather than duplicating full new
references.

## 7. Developer architecture documentation added

`docs/dev/synthesis_engine.md` records the production lifecycle, cache
invalidation, exact batched gate, interval-index mechanism, immutable line
state, plane-parallel/spherical split, ray-order correction, irregular-grid
resampling, numerical status, adopted memory architecture, dynamic transfer
capacity, and process-global concurrency model.

Rejected prototypes are not described as architecture.

## 8. Release notes and changelog

`analysis/v120_release_notes_draft.md` is a non-published GitHub release-note
draft with highlights, migration table, non-multiplicative component results,
correctness changes, numerical validation, memory, compatibility, and current
limitations.

The root changelog now separates performance, correctness, and compatibility,
and uses the restrained 19--40x canonical 10 A headline. The documentation
changelog adds only a concise v1.2 release-candidate summary and links to the
migration page.

## 9. Unresolved documentation questions

1. The release is still untagged, so the root changelog remains under
   `Unreleased` and the GitHub release text remains an analysis draft.
2. Legacy controls cannot reproduce every v1.1 result because the irregular-
   grid integration and spherical ray-order corrections remain active. The
   docs state this rather than promising bitwise reproduction.
3. `SME_DLL.session()` exists but has no established high-level public-API
   commitment. It is intentionally omitted from user guidance.
4. No combined 800 A end-to-end RSS run was repeated after the streaming parser
   change. Documentation quotes the separately measured parser and native-cache
   reductions, not an invented combined peak.
5. Nineteen pre-existing Sphinx warnings remain outside this task: old title
   underlines, a `concepts/lsf` versus `lfs` toctree typo, unlisted legacy docs,
   and old ambiguous/missing MyST references. None originates in the new or
   edited v1.2 links.
6. `line_select_almax_threshold` is not a membership-only override: source
   passes it into `ALMAXRange`, so it also affects the ranges produced by that
   precomputation. The documentation now records this coupling rather than
   claiming that it cleanly separates selection from support.

## 10. Documentation build result

Command:

```text
/Users/mingjie/opt/miniconda3/envs/astro_sme/bin/python -m sphinx \
    -b html docs /tmp/pysme-v120-docs-html
```

Result: HTML build succeeded for 51 source files. Sphinx reported 19 existing
warnings listed above; the new migration, performance, architecture, changelog,
and cross-page links introduced no warning.

The release-candidate visual preview was rebuilt with
`-D release=1.2.0 -D version=1.2`, so its homepage shows `Version: 1.2.0`
without hard-coding a pre-release version in `docs/conf.py`. This full preview
build also succeeded with the same 19 existing warnings.

## 11. Release-facing claims and benchmark sources

| Claim | Source |
|---|---|
| Canonical 10 A complete synthesis improved about 19--40x | `analysis/v120_release_performance_summary.md`, Tables A/B |
| Tested Mg b/H-alpha complete synthesis improved about 18--43x | `analysis/v120_release_performance_summary.md`, Table A |
| Continuum/CDR task improved about 17.2x; hundreds-fold exact-query reduction | `analysis/v120_release_performance_summary.md`, sections 6--7 |
| Direct-synthesis continuum effect at most about `1.4e-9` | `analysis/v120_release_performance_summary.md`, section 7.1 |
| Cumulative-bin CDR boundary case about `7.24e-6` | same source, section 7.1 |
| Batching controlled outputs bitwise identical | same source section 7.2; plane-parallel and spherical batched audits |
| Candidate visits reduced about 99.63--99.94% | same source, section 6 |
| Float cache reduced three-array payload 1.66 GB to 0.83 GB and synthesis RSS by about 0.81 GB | `analysis/v120_line_state_memory_audit.md` |
| Float cache maximum tested flux difference about `1.6e-8`, unchanged masks/ranges | same source |
| Streaming parser peak decreased about 56% for 1.26 million lines | `analysis/v120_vald_parser_memory_audit.md` |
| Release-candidate suite: 187 passed, 1 expected failure | `analysis/v120_release_performance_summary.md` and line-state audit |

## Statements intentionally not documented as v1.2 features

- rejected `VVOIGT` factorization;
- wavelength-chunked line-state materialization;
- ALMAX-seeded adaptive-grid and fixed regular velocity-grid prototypes;
- TiO molecular cross sections;
- any broader “PreparedSynthesis” abstraction beyond the implemented
  abundance-only state reuse;
- true thread-parallel or instance-reentrant SMElib;
- internal-selection batching;
- a hard claim for combined post-parser end-to-end peak memory;
- `long_continuum` as an `SME_Structure` setting;
- `SME_DLL.session()` as a supported general user parallelism API.

## Page summary

| File/page | Audience | Main change | Release-specific or permanent |
|---|---|---|---|
| `docs/getting_started/whats_new_v120.md` | All users | v1.2 overview, result expectations, and essential migration guidance | Release-specific |
| `docs/index.rst` | All users | v1.2 highlight and overview link | Release-specific highlight |
| `docs/advance/synthesis_performance.md` | Users/advanced users | Current path table, algorithm concepts, numerical controls | Permanent |
| `docs/advance/line_filtering.md` | Advanced users | Membership/support/sampling separation | Permanent |
| `docs/advance/line_selection_reference.md` | Advanced users | Exact current semantics and `wave`/`wint` distinction | Permanent |
| `docs/fundamentals/flux_inten.md` | Users/advanced users | Native irregular grid and corrected integration order | Permanent |
| `docs/concepts/sme_struct.md` | Users | Correct public field meanings/defaults | Permanent |
| `docs/fundamentals/linelist.md` | Users | Incremental VALD parsing and memory note | Permanent |
| `docs/advance/faq.md` | Users | Supported parallelism guidance | Permanent |
| `docs/dev/synthesis_engine.md` | Developers | Production architecture, validation, memory, concurrency | Permanent |
| `docs/dev/changelog.md` and `changelog.md` | Upgraders/maintainers | Concise v1.2 summary | Release-specific |
| `analysis/v120_release_notes_draft.md` | Release maintainers | Draft GitHub release notes | Release-specific |
| `README.md` | Prospective users | One concise adaptive-synthesis feature statement | Permanent |
