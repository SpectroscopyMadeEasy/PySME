# PySME v1.2

**Faster synthesis, lower memory use, and improved numerical robustness**

PySME v1.2 makes standard spectral synthesis substantially faster.
Representative calculations in our validation set are typically **20–40×
faster**, and most existing synthesis and fitting scripts benefit automatically
without any code changes.

The release also reduces memory use for very large line lists and includes
several numerical-correctness improvements. The gains come from avoiding
redundant calculations while retaining the same radiative-transfer treatment
and accuracy controls.

## What is new?

PySME v1.2 removes major sources of unnecessary computation from spectral
synthesis:

- continuum-opacity calculations are cached and reused instead of being
  repeated for many lines and wavelengths;
- adaptive synthesis evaluates only the spectral lines that can contribute at
  each wavelength, rather than repeatedly checking the full line list.

The improvement is particularly large for line-rich spectra.

The release also reduces memory use for very large VALD line lists and fixes
several long-standing numerical edge cases.

## Do I need to do anything?

**Usually, no.** Existing synthesis and fitting scripts generally benefit
automatically after upgrading to v1.2.

There is no general requirement to rerun LTE analyses solely because of the
v1.2 performance improvements.

## Will my scientific results change?

For most calculations, differences are very small. The performance
improvements do not affect radiative-transfer accuracy.

Some spectra can differ slightly from older PySME versions because v1.2 also
corrects several historical numerical behaviours. In our validation set, most
differences were at the level of about 10⁻⁵ or smaller in normalized flux;
the largest case was about 10⁻³, in a metal-poor spectrum affected by the
historical weak-line pruning. Where the differences were largest, the v1.2
results agreed better with reference syntheses in which substantially more
weak lines were retained.

If exact comparison with older PySME output is important, compatibility
controls are available for reproducing the historical synthesis paths where
supported. See {ref}`compatibility and numerical controls
<compatibility-controls>`.

See {ref}`numerical changes and validation <numerical-behavior-and-validation>`
for details.

## Find the level of detail you need

- **[Performance and numerical controls](../advance/synthesis_performance.md)** —
  synthesis paths, validation, compatibility options, and accuracy controls.
- **[Release notes](../dev/changelog.md)** — the complete high-level list of
  changes.
