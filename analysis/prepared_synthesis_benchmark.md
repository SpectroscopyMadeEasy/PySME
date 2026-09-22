# Prepared synthesis state-reuse benchmark

## Question

Measure the ceiling for repeated abundance synthesis when an unchanged line
list and atmosphere/model remain resident in SMElib.  The resulting
optimization is used automatically inside abundance-only `solve(...)`
lifecycles; no public `PreparedSynthesis` object is added.

The reference path repeats the current sequence on every evaluation:

```text
InputLineList -> InputModel -> InputAbund -> Ionization -> Opacity -> Transf
```

The prepared path performs the full sequence once, then repeats only:

```text
InputAbund -> Ionization -> Opacity -> Transf
```

Both paths enable the existing exact EOS warm start.  The prepared path does
not reuse EOS output, continuum opacity, or line opacity across abundance
changes.

## Inputs

- Line list: `merged_3700-9500_hfs.lin`, extracted with about 150 A padding.
- Centre: 5200 A.
- Output sampling: 0.05 A.
- Windows: 10 A (378,154 input lines) and 200 A (592,757 input lines).
- Models: solar dwarf (5777 K, log(g)=4.44, [M/H]=0.0) and cool metal-rich
  dwarf (4250 K, log(g)=4.5, [M/H]=+0.3).
- Six Fe evaluations: offsets 0.00, +0.05, -0.05, +0.10, -0.10, 0.00 dex.
- `accrt=1e-4`, ALMAX selection, adaptive continuum grid.
- One full preparation call is excluded from each arm's evaluation timings.

## Results

Times are medians across the six abundance evaluations.  Total speedup is the
ratio of the sums of all six wall times.

| model | window | current | prepared | median speedup | total speedup | max abs flux delta |
|---|---:|---:|---:|---:|---:|---:|
| solar dwarf | 10 A | 0.561 s | 0.373 s | 1.50x | 1.52x | 0 |
| solar dwarf | 200 A | 0.956 s | 0.684 s | 1.40x | 1.38x | 0 |
| cool metal-rich dwarf | 10 A | 2.018 s | 1.833 s | 1.10x | 1.09x | 0 |
| cool metal-rich dwarf | 200 A | 4.425 s | 4.144 s | 1.07x | 1.07x | 0 |

## Interpretation

State reuse removes a roughly 0.19--0.28 s fixed cost per evaluation in these
cases.  That is significant for the solar model, where it reduces repeated
synthesis wall time by about 29--34%.  It is only a 6--9% reduction for the
cool metal-rich model because line opacity and transfer dominate its runtime.

The spectra were bit-for-bit identical to the current repeated-input path for
all 24 A/B comparisons.  This is expected: the prototype changes lifecycle and
invalidation only, not physics.  It reuses the same ALMAX information that the
current `line_select_recompute="if_stale"` path already reuses for abundance
changes; it is therefore an equivalence test against current behaviour, not an
independent validation of ALMAX staleness for arbitrarily large abundance
changes.

The result supports the implemented, narrowly scoped solver optimization for
abundance-only iterations.  It does not support treating `PreparedSynthesis`
as a large general-purpose state API for cool, line-rich spectra.  PySME keeps
the full setup path for mixed atmosphere/abundance fits, dynamic line-list
selection, forced line-info recomputation, callable atmospheres, and line-list
parameter fits.

## Artifacts

- Driver: `analysis/prepared_synthesis_benchmark.py`
- Raw 10 A output: `analysis/prepared_synthesis_benchmark_10A.json`
- Raw 200 A output: `analysis/prepared_synthesis_benchmark_200A.json`
