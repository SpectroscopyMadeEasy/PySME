# Prepared synthesis state-reuse benchmark

## Question

Measure the ceiling for repeated abundance synthesis when an unchanged line
list and atmosphere/model remain resident in SMElib.  The resulting
optimization is used automatically inside abundance-only `solve(...)`
lifecycles; no public `PreparedSynthesis` object is added.

The reference path repeats the current sequence on every evaluation:

```text
InputLineList -> InputModel -> InputAbund -> Ionization -> Opacity
              -> ALMAXRange -> InputLinePrecomputedInfo -> Transf
```

The prepared path performs the full sequence once, then repeats only:

```text
InputAbund -> Ionization -> Opacity
           -> ALMAXRange -> InputLinePrecomputedInfo -> Transf
```

Both paths enable the existing exact EOS warm start.  The prepared path does
not reuse ALMAX/range metadata, EOS output, continuum opacity, or line opacity
across abundance changes.  `Transf` can consume the opacity/Voigt state from
the immediately preceding `ALMAXRange`, avoiding a duplicate `LINEOPAC` pass.

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
| solar dwarf | 10 A | 2.046 s | 2.094 s | 0.98x | 1.00x | 0 |
| solar dwarf | 200 A | 3.132 s | 2.903 s | 1.08x | 1.12x | 0 |
| cool metal-rich dwarf | 10 A | 2.854 s | 2.595 s | 1.10x | 1.09x | 0 |
| cool metal-rich dwarf | 200 A | 5.680 s | 5.404 s | 1.05x | 1.08x | 0 |

## Interpretation

Once ALMAX is correctly invalidated by abundance changes, state reuse is a
small setup optimization: 0--12% over these complete six-evaluation workloads.
At 10 A the solar difference is within run-to-run noise; at 200 A it saves
about 8--12%.  The cool model saves about 5--10% because line opacity and
transfer dominate its runtime.

The spectra were bit-for-bit identical to the current repeated-input path for
all 24 A/B comparisons.  This is expected: the optimization changes lifecycle,
not physics, and both arms recompute abundance-dependent ALMAX information.

The result supports keeping the implemented, narrowly scoped solver
optimization, but it is not a major speed feature after correct ALMAX
invalidation.  PySME keeps the full setup path for CDR selection, mixed
atmosphere/abundance fits, dynamic line-list selection, forced line-info
recomputation, callable atmospheres, and line-list parameter fits.

## Artifacts

- Driver: `analysis/prepared_synthesis_benchmark.py`
- Raw 10 A output: `analysis/prepared_synthesis_benchmark_10A.json`
- Raw 200 A output: `analysis/prepared_synthesis_benchmark_200A.json`
