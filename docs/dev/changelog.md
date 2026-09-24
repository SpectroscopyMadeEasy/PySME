# Changelog

This page keeps only high-level release and migration notes.

Detailed per-release changes are maintained on GitHub Releases:

- https://github.com/SpectroscopyMadeEasy/PySME/releases
- https://github.com/SpectroscopyMadeEasy/PySME/tags

## Changelog policy

- Use GitHub Releases as the source of truth for detailed change history (v0.7.x and newer).
- Use the Git tags list for older historical versions.
- Keep this page concise, focused on major milestones and migration notes.

## v1.2.0 (release candidate)

- Reworked adaptive synthesis with continuum-opacity caching and
  generation-batched, interval-indexed transfer for plane-parallel and
  spherical atmospheres.
- Reduced memory use for large line lists through float line-state caches and
  incremental parsing of supported large VALD files.
- Corrected historical weak-line pruning, mutable line-range state,
  irregular-grid flux integration, spherical ray ordering, and fixed
  adaptive-transfer capacity behavior.
- Complete synthesis was approximately 19--40x faster in the canonical 10 A
  validation matrix; gains depend on the workload.

Existing scripts generally require no changes. See [](../getting_started/whats_new_v120.md)
for the user overview and [](../advance/synthesis_performance.md) for
optimized-path conditions, numerical differences, and compatibility controls.

## v1.1.1

- Fixed a critical NLTE line-indexing error when SMElib discards transitions
  with unsupported ionization stages. PySME now maps Python line-list indices
  to the compact internal SMElib indices before assigning departure
  coefficients.
- The affected releases are v0.4.151 through v1.1.0. These versions are no
  longer recommended for scientific NLTE synthesis, and affected calculations
  should be rerun with v1.1.1 or later.
- This is a Python-side correctness fix and does not require a new SMElib
  release.

## v1.1.0

- Added optional continuum scattering, strict NLTE fallback handling, improved spherical-atmosphere interpolation, and updated solar-abundance and data-download support.
- PySME v1.1.0 is paired with SMElib v6.13.19.

## Legacy 0.x notes

Historical details for older 0.x releases are available in the tag list.
