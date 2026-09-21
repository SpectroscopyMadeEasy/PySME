# Production adaptive continuum-grid benchmark

This follow-up tests the production SMElib implementation, rather than the
temporary fixed-grid instrumentation used in
`continuum_opacity_grid_benchmark.md`.

## Configuration

- nominal spacing: 1 A
- tolerances: `3e-3`, `1e-3`, `3e-4`, `1e-4`
- minimum recursive spacing: `1e-3 A`
- exact guard band: `+/-0.02 A` around known H I, Mg I, and Si I edges
- cached state: all 13 continuum opacity-source components
- CDR line list: 378,154 valid lines over 5045--5355 A
- models: solar dwarf, metal-poor dwarf, cool metal-rich dwarf, and giant

Machine-readable results are in `adaptive_continuum_grid_benchmark.json`.

## Continuum accuracy over 3500--9500 A

For the solar model, using 16,040 exact reference wavelengths:

| tolerance | exact nodes/evaluations | refined intervals | p99 relative error | max relative error |
|---:|---:|---:|---:|---:|
| `3e-3` | 24,511 / 24,590 | 2 | `6.77e-6` | `1.27e-3` |
| `1e-3` | 24,536 / 24,615 | 9 | `6.16e-6` | `5.03e-4` |
| `3e-4` | 24,607 / 24,686 | 30 | `3.52e-6` | `2.67e-4` |
| `1e-4` | 24,705 / 24,784 | 62 | `1.42e-6` | `9.27e-5` |

The high node count in this diagnostic scan is expected: its 0.5 A sampling
queries the quarter/midpoint certification nodes across the full 6000 A
range. The CDR workload below only materializes 1,246--1,420 nodes.

## CDR runtime and selection

At the default `1e-3` tolerance, ALMAX/range plus central-depth computation
took 3.56 s (solar), 3.62 s (metal-poor), 4.37 s (cool metal-rich), and
4.07 s (giant). A new exact solar run using the identical `ALMAXRange +
CentralDepth` call chain took 61.23 s, so the directly comparable solar speedup
is 17.2x. The earlier 90.8--100.6 s references included a somewhat different
`Transf + CentralDepth` workflow and should not be divided directly by these
new times. The production grid uses 1,246--1,380 exact continuum nodes for
0.76--1.14 million continuum queries.

Changing the tolerance did not change the final cumulative-bin line-selection
differences relative to the exact reference:

| model | false positive | false negative | max normalized-flux difference |
|---|---:|---:|---:|
| solar dwarf | 1 | 1 | `7.13e-10` |
| metal-poor dwarf | 1 | 2 | `6.80e-10` |
| cool metal-rich dwarf | 16 | 12 | `7.24e-6` |
| giant | 1 | 0 | `5.71e-10` |

The cool-model difference is the previously identified cumulative selection
boundary effect, not a direct continuum interpolation error. Tightening the
continuum tolerance does not change those same weak-line rank decisions.

## Mg I 3756.607779 A regression

An initial recursive implementation with `min_step=1e-3 A` still left a
localized `~2e-5` normalized-flux error immediately redward of the Mg I edge.
Adding an exact `+/-0.02 A` guard band reduced the 3756.55--3756.70 A test to:

- max `|delta F_norm| = 3.58e-11`
- RMS `= 7.47e-12`

This result is independent of the four tested tolerances because the dangerous
part of the edge shoulder is deliberately evaluated exactly.

## Decision

`rtol=1e-3` is the recommended default. It keeps the continuum error below the
requested scale in the broad scan, gives the same line-selection and spectrum
results as tighter tolerances in the tested CDR workload, and avoids spending
extra work where it has no measurable spectrum benefit. The exact mode remains
available for reference calculations.
