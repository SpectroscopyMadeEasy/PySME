# PySME v1.2.0 synthesis-performance release audit

## 1. Executive summary

**Recommendation: GO for a v1.2.0 release candidate.**

The two production optimizations are independently validated and materially
improve complete `synthesize_spectrum()` wall time.  On the six canonical
10 A workloads, the historical-to-new end-to-end speed-up is **19.4--39.8x**.
The 25 A Mg b and H-alpha workloads improve by **17.9--43.2x**.  New-path
200 A synthesis completes in 4.0 s for the Sun and about 46--50 s for the two
line-rich giant/cool cases; 800 A completes in 11.3 s for the Sun, 8.24 s for
the metal-poor giant, and 161--171 s for line-rich giant/cool cases.

The controlled 2x2 ablation shows that the gains are not multiplicative
component headlines.  Continuum caching alone gives 11.9--17.2x for ordinary
solar/metal-poor cases but only about 1.5x where adaptive transfer dominates.
Batching alone gives about 2.8--3.0x end-to-end for line-rich cases but little
benefit when exact continuum preparation dominates.  Together they give
19.3--36.7x against the semantics-controlled baseline.  The measured positive
interaction is 1.1--8.2x relative to the naive product.

There is no batching approximation under controlled semantics: C0B0 and C0B1
post-integration spectra are bitwise identical in all six end-to-end 10 A
cases, and the dedicated native audits show bitwise-identical common-node
`sint`/`cint`, including spherical grazing rays and an H-NLTE departure-
coefficient case.  Adaptive continuum interpolation changes normalized flux
by at most `1.40e-9` in this direct synthesis matrix.  The separate cumulative-
bin CDR audit has a larger worst case, `7.24e-6`, caused by weak-line rank
decisions at the selection boundary rather than the transfer calculation.

The audit found and fixed one release blocker: a cool, line-rich 800 A segment
could exceed the historical 400,000-point transfer allocation.  High-level
synthesis now estimates capacity from the in-window line-centre upper bound
and the RKINTS 0.3 km/s spacing floor.  Ordinary windows retain the 400,000
default; the formerly failing 800 A case now completes in 170.69 s.  The full
suite passes: **187 passed, 1 xfailed, 15 warnings**.

## 2. What changed

The release combines performance and correctness work that must not be
described as one numerical approximation:

| Change | Purpose | Numerical/physical effect |
|---|---|---|
| Adaptive continuum-opacity cache | Replace repeated exact `CONTOP` calls with an edge-aware adaptive cache | Interpolation approximation; max direct-synthesis `|dFnorm|=1.40e-9` here |
| Generation-batched adaptive RKINTS | Expose sorted generations to `LineIntervalSweep + OPMTRXIndexed` | Scheduling only; controlled native intensities are bitwise identical |
| Remove sequential `MARK=2` second pruning | Keep ALMAX membership immutable | Intended semantic/correctness change; closer to minimally pruned reference |
| Immutable physical `Wlim` | Prevent synthesis from shortening precomputed support | Intended correctness change; removes grid-layout/state-history dependence |
| Irregular-grid flux integration fix | Resample every mu intensity to a common regular log-lambda grid before disk integration/broadening | Correctness fix; removes false neighbour/grid dependence |
| Spherical ray-order fix | Evaluate normal rays before deepest-to-shallowest grazing rays | Correctness fix; removes caller mu-order dependence |
| Adaptive transfer capacity sizing | Avoid a fixed 400k ceiling for unusually wide line-rich segments | Capacity/robustness only; no physics change |

The line semantics are now explicit:

```text
ALMAX                         -> line membership
accrt + physical Wlim         -> wavelength support
accwi                         -> adaptive transfer-grid sampling accuracy
```

Production batched transfer does not apply the historical sequential
`MARK=2` pruning and does not mutate precomputed physical ranges.  Explicit
legacy mode and internal line selection retain legacy RKINTS.

## 3. Benchmark environment and hygiene

- Machine: Apple M5, 10 logical cores, macOS 26.6.2 arm64.
- Python 3.11.14, NumPy 2.4.6, single process, no line-selection workers,
  `OMP_NUM_THREADS` unset.
- GCC/GFortran 16.1.0; release library built with `g++-16 -O3` plus the
  configured `-g -O2` flags.
- Parent commit: `8b6a994094abbba73a161bac123de181d04e0201`.
- SMElib commit: `c92b466be46a94abf97fb39c19dd948af149f149`.
- Measured working-tree binary-diff hashes are stored in
  `v120_release_benchmark.json`; the parent and SMElib values were
  `c8e2d46f...35c3` and `4818c34a...2a32`.  Later report-only edits
  do not affect the measured binaries.
- `accrt=1e-4`, `accwi=3e-3`; continuum base step 1 A, `rtol=1e-3`, minimum
  step 0.001 A; output grid 0.05 A; `R=60,000`, zero rotation/macroturbulence.
- VALD parsing is excluded.  Every timing uses a fresh `SME_Structure` and
  `Synthesizer`; the OS file cache is warm.
- Fast cases use at least three runs and report medians.  Historical and
  controlled line-rich cases taking several minutes use one run, as marked.
- Main line-list SHA256 values: 10 A `03c90a90...993e`, 200 A
  `02ff0142...49b9`, 800 A `e79601e1...de09`, H-alpha
  `8de1d77f...3831`.
- Spherical atmosphere SHA256 values: moderate giant `54928d6d...3ecc`, cool
  giant `7586c5df...d8e1`, metal-poor giant `0e6e6539...9249`.

Historical H means exact continuum, internal line selection, legacy
RKINTS/RKINTS_sph, and historical mutable MARK/Wlim behavior.  Controlled C0B0
uses exact continuum, an immutable precomputed ALMAX mask/ranges, no second
pruning, and sequential legacy scheduling.  A temporary non-public audit patch
constructed that benchmark baseline.  The patch and its environment switch
were removed after the controlled measurements, as required for release
hygiene; they never affected H or production C1B1 when unset.

## 4. Table A: release/user-facing complete synthesis

| Workload | Geometry | Historical H (s) | New C1B1 (s) | Speed-up | Max `|H-new|` |
|---|---|---:|---:|---:|---:|
| Solar dwarf, 10 A | PP | 41.141 (n=3) | 1.902 (n=4) | **21.63x** | `1.83e-5` |
| Metal-poor dwarf, 10 A | PP | 35.584 (n=3) | 1.774 (n=3) | **20.06x** | `1.34e-3` |
| Cool metal-rich dwarf, 10 A | PP | 238.685 (n=1) | 5.995 (n=4) | **39.82x** | `1.40e-9` |
| Moderate giant, 10 A | spherical | 187.570 (n=1) | 5.600 (n=4) | **33.50x** | `4.46e-8` |
| Cool line-rich giant, 10 A | spherical | 212.600 (n=1) | 6.089 (n=3) | **34.92x** | `9.69e-10` |
| Metal-poor giant, 10 A | spherical | 36.359 (n=3) | 1.872 (n=3) | **19.42x** | `4.88e-6` |
| Solar Mg b, 25 A | PP | 40.236 (n=1) | 1.917 (n=3) | **20.99x** | `1.06e-5` |
| Moderate-giant Mg b, 25 A | spherical | 298.986 (n=1) | 8.557 (n=3) | **34.94x** | `3.37e-8` |
| Solar H-alpha, 25 A | PP | 52.610 (n=1) | 2.934 (n=3) | **17.93x** | `3.42e-6` |
| Moderate-giant H-alpha, 25 A | spherical | 422.132 (n=1) | 9.780 (n=3) | **43.16x** | `5.41e-8` |

The H-to-new differences in the last column are not batching error.  Depending
on the model they include removal of historical second pruning, immutable
ranges, irregular-flux integration, the spherical ray-order correction, and
continuum interpolation.  The controlled comparisons below isolate the two
performance optimizations.

For wider production-only checks (historical wide cases were intentionally not
run because of their multi-minute/hour cost):

| Workload | Geometry | New production median |
|---|---|---:|
| Solar, 200 A | PP | 3.984 s (n=4) |
| Cool metal-rich dwarf, 200 A | PP | 50.475 s (n=4) |
| Moderate giant, 200 A | spherical | 46.442 s (n=4) |
| Solar, 800 A | PP | 11.251 s (n=4) |
| Metal-poor giant, 800 A | spherical | 8.239 s (n=3) |
| Cool metal-rich dwarf, 800 A | PP | 170.685 s (n=1) |
| Moderate giant, 800 A | spherical | 160.801 s (n=1) |

These runs show that the production path scales beyond the 10 A microbenchmark.
The older isolated v1.1.1 comparison already established a `>104x` lower
bound for a 200 A fixed-grid solar synthesis stopped after 360 s; it is not
mixed into the adaptive-grid table above.

## 5. Table B: controlled 2x2 attribution

| 10 A workload | C0B0 (s) | C1B0 (s) | C0B1 (s) | C1B1 (s) | Continuum only | Batch only | Combined | Interaction |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Solar dwarf | 37.615 | 3.170 | 35.705 | 1.902 | 11.86x | 1.05x | **19.78x** | 1.58x |
| Metal-poor dwarf | 35.690 | 2.081 | 34.540 | 1.774 | 17.15x | 1.03x | **20.12x** | 1.14x |
| Cool metal-rich dwarf | 210.437 | 140.153 | 74.186 | 5.995 | 1.50x | 2.84x | **35.10x** | 8.23x |
| Moderate giant | 187.843 | 120.741 | 66.740 | 5.600 | 1.56x | 2.81x | **33.54x** | 7.66x |
| Cool line-rich giant | 223.335 | 149.356 | 73.313 | 6.089 | 1.50x | 3.05x | **36.68x** | 8.05x |
| Metal-poor giant | 36.052 | 2.339 | 36.298 | 1.872 | 15.42x | 0.99x | **19.26x** | 1.26x |

`Interaction = S_total / (S_cont * S_batch)`.  Exact continuum preparation
dominates C0B1 in solar/metal-poor cases, while a sequential full-scan transfer
dominates C1B0 in line-rich cases.  Removing both bottlenecks therefore exposes
a much smaller residual than either one-factor run, producing the strong
positive interaction.  It would be wrong to quote the isolated factors as a
simple product.

## 6. Stage and counter attribution

Representative first-run stage timings make the bottlenecks explicit:

- Solar C0B0 spends about 35.0 s in ALMAX/range preparation and 1.82 s in
  transfer.  C1B1 reduces these to about 1.35 s and 0.09 s.
- Cool-dwarf C0B0 spends 66.1 s in ALMAX/range preparation and 143.5 s in
  transfer.  C1B1 reduces these to about 1.68 s and 3.72 s.
- Moderate-giant C0B0 spends 64.7 s in ALMAX/range preparation and 122.5 s in
  transfer.  C1B1 reduces these to about 1.68 s and 3.14 s.
- Atmosphere interpolation, EOS, flux integration, and broadening are included
  in all complete times.  They are not the dominant cost in this matrix.

Continuum counters for C1B1 include:

| Workload | Continuum queries | Exact evaluations/nodes | Query-to-exact reduction |
|---|---:|---:|---:|
| Solar 10 A | 399,538 | 1,254 | 318.6x |
| Cool dwarf 10 A | 775,509 | 1,380 | 562.0x |
| Moderate giant 10 A | 723,131 | 1,332 | 542.9x |
| Cool dwarf 800 A | 3,077,247 | 4,745 | 648.5x |

The canonical standalone continuum/CDR benchmark is more directly comparable:
exact `ALMAXRange + CentralDepth` took 61.23 s for the Sun, while adaptive
continuum took 3.56 s, a **17.2x** speed-up.  The four CDR models reduced
0.76--1.14 million continuum queries to about 1,246--1,380 exact nodes.

Batching counters show that it does not reduce physical transfer evaluations.
For the 10 A cool dwarf, both schedulers evaluate 10,049 wavelengths; the
indexed path visits 8.59 million candidates instead of a 2.524 billion
full-scan equivalent, a **99.66% reduction**.  Across this release matrix the
reduction is 99.63--99.94%.  Native `TBINTG`/`TBINTG_sph` time is nearly
unchanged; `OPMTRX` candidate traversal is the removed bottleneck.

Dedicated native headline results are:

- Plane-parallel cool metal-rich transfer: 115.23 s legacy versus 4.82 s
  production, **23.91x**.
- Spherical dense cases: **31.35x**, **23.88x**, and **6.88x** for moderate,
  cool, and metal-poor giants.
- Spherical Mg b and H-alpha: **39.46x** and **49.83x**.

## 7. Accuracy audit

### 7.1 Adaptive continuum cache: C0B0 versus C1B0

| 10 A workload | Max `|dFnorm|` | RMS | p99 | `|dEW|` (A) | Mask differences | Range differences |
|---|---:|---:|---:|---:|---:|---:|
| Solar dwarf | `7.83e-10` | `2.36e-10` | `7.68e-10` | `1.33e-9` | 0 | 0 |
| Metal-poor dwarf | `7.47e-10` | `1.16e-10` | `4.23e-10` | `2.82e-10` | 0 | 0 |
| Cool metal-rich dwarf | `1.40e-9` | `6.63e-10` | `1.31e-9` | `5.86e-9` | 0 | 0 |
| Moderate giant | `5.44e-10` | `1.12e-10` | `3.25e-10` | `7.15e-10` | 0 | 0 |
| Cool line-rich giant | `9.69e-10` | `5.12e-10` | `9.13e-10` | `4.53e-9` | 0 | 0 |
| Metal-poor giant | `4.85e-10` | `8.30e-11` | `4.11e-10` | `2.32e-10` | 0 | 0 |

Raw ALMAX ratios are not bitwise equal because their denominator uses the
interpolated continuum; the largest absolute differences range from
`2.42e-5` to 0.490 and occur for already-strong ratios.  The decision-relevant
membership masks and all physical range endpoints are identical in this
matrix.

The broad continuum validation over 3500--9500 A gave p99 relative continuum
error `6.16e-6` and maximum `5.03e-4` at the default `rtol=1e-3`.  Explicit
physical knots and the exact +/-0.02 A edge guard reduce the Mg I 3756.608 A
regression to max `|dFnorm|=3.58e-11`, RMS `7.47e-12`.  H I thresholds are
also explicit knots/guarded regions.

The separate cumulative-bin CDR selection audit found 1/1, 1/2, 16/12, and
1/0 false-positive/false-negative changes for solar, metal-poor, cool, and
giant cases, with maximum flux differences `7.13e-10`, `6.80e-10`,
`7.24e-6`, and `5.71e-10`.  The cool result is a cumulative weak-line ranking
boundary; tightening continuum tolerance did not change those decisions.

### 7.2 Generation batching: C0B0 versus C0B1

All six complete 10 A spectra are bitwise equal, with zero flux and EW
difference.  More stringent native checks show:

- accepted common-node wavelength arrays and every `sint(mu,lambda)` and
  `cint(mu,lambda)` are bitwise identical to fixed-grid evaluation;
- the immutable ALMAX mask and every physical range remain bitwise unchanged;
- spherical equality includes the two grazing rays, intermediate rays, and
  disk-centre ray;
- the validated H-NLTE departure-coefficient case is bitwise identical under
  controlled semantics;
- hydrogen special profiles remain inside the unchanged `OPMTRX` physics path.

Thus batching changes scheduling and candidate lookup only.  The existing
`accwi=3e-3` adaptive-grid interpolation criterion remains a separate legacy
numerical heuristic and is not a global spectrum-error guarantee.

### 7.3 Combined and historical comparisons

C0B0 versus C1B1 has exactly the same flux differences as C0B0 versus C1B0,
because controlled batching is exact.  Historical H versus C1B1 can be much
larger: `1.34e-3` for the metal-poor dwarf, for example.  The second-pruning
audit demonstrated that this is the removed order-dependent `MARK=2` behavior,
which moved the result away from a minimally pruned reference.  Solar/giant
`~1e-5--1e-7` historical residuals were traced to mutable Wlim state.

These are intended correctness/semantic changes, not performance-approximation
errors.  The spherical ray-order and irregular-flux fixes are likewise
correctness changes and must be described separately in release notes.

## 8. H/NLTE and special profiles

Both plane-parallel and spherical H-alpha complete-synthesis cases are included
in Table A.  The dedicated spherical H-NLTE audit uses `nlte_H_pysme.grd` and
keeps the existing departure-coefficient opacity/source calculation inside
`OPMTRX`; controlled common-node results are bitwise identical.  Its smaller
1.35x performance gain is expected for the compact 983-line regression list
and is not used as the line-rich headline.

Mg b is represented in both geometries.  The end-to-end spherical case improves
34.94x, while the dedicated native transfer audit improves 39.46x.  No special
line-profile branch was duplicated or bypassed.

## 9. Solve timing

The current solar, wave-only, one-parameter Fe solve was rerun after batching:

- 10 synthesis calls, five residual iterations;
- total solve 17.114 s, of which 17.094 s is synthesis;
- 1.709 s per synthesis call;
- recovered Fe to `1.22e-7 dex` and final max flux residual `2.13e-7`.

The corresponding v1.1.1 solve had previously been stopped after 1,011 s
without completing, so the current result establishes a conservative **>59x**
end-to-end lower bound.  The completed fixed-grid solve comparison remains
26.69x solar, 24.13x metal-poor, and 14.28x cool metal-rich; those figures also
include EOS warm-start and prepared abundance-state reuse and are not isolated
batching results.

## 10. Memory sanity

Each row below is a fresh-process production measurement.  Peak RSS includes
the parsed VALD/Pandas line list, its synthesis copy, atmospheres, Python, and
the native library; it is not the cache size.

| Workload | Peak RSS | Continuum payload | Generation payload | Transfer nodes |
|---|---:|---:|---:|---:|
| Solar 10 A | 2.17 GiB | 6.8 MiB | 0.3 MiB | 2,447 |
| Cool dwarf 10 A | 2.19 GiB | 7.5 MiB | 1.2 MiB | 10,049 |
| Moderate giant 10 A | 2.18 GiB | 7.4 MiB | 1.1 MiB | 9,885 |
| Solar 200 A | 3.19 GiB | 11.4 MiB | 4.0 MiB | 34,561 |
| Cool dwarf 200 A | 3.34 GiB | 11.8 MiB | 15.7 MiB | 137,157 |
| Moderate giant 200 A | 3.22 GiB | 11.9 MiB | 15.4 MiB | 134,933 |
| Solar 800 A | 3.64 GiB | 24.8 MiB | 14.7 MiB | 128,379 |
| Cool dwarf 800 A | 4.45 GiB | 25.0 MiB | 60.9 MiB | 532,481 |
| Moderate giant 800 A | 4.60 GiB | 25.3 MiB | 59.3 MiB | 518,361 |

The 10-to-800 A RSS growth is mainly the 378k-to-1.26M-line input and Python
copies, not a continuum-cache explosion.  Algorithm-specific payload remains
about 85 MiB in the worst measured case.  The dynamic capacity allocation adds
storage only when the estimated adaptive grid can exceed 400k.  This is
healthy for the requested workloads, though single-segment line-rich 800 A
synthesis is intentionally documented as a high-memory use case.

## 11. Continuum-cache invalidation

Source audit and new regression coverage establish:

- `InputModel` clears the cache, covering atmosphere, geometry, and model
  opacity flags;
- `InputAbund` clears it;
- `Ionization` clears it, covering EOS/species/electron-state changes;
- changing continuum-grid mode/tolerance/spacing clears it;
- line-list-only input preserves it because line data do not enter `CONTOP`;
  the normal following `Ionization` still clears it when the species/EOS state
  can change.

The optional continuum-scattering source mode consumes the same cached 13
physical opacity components, so changing source treatment does not require
recomputing those components.  H2 broadening is a line-broadening setting and
does not alter continuum opacity.

## 12. Known fallback paths

The production batched path is used for plane-parallel and spherical adaptive
transfer when long continuum and validated precomputed ALMAX/CDR mask/ranges
are present.  The following intentionally retain established paths:

- user-supplied fixed wavelength grids (already evaluated as a sorted grid);
- `line_select_method="internal"`;
- `long_continuum=False`;
- explicit `transfer_grid_method="legacy"`;
- missing or stale precomputed line information.

No PreparedSynthesis expansion, TiO cross-section work, threaded/re-entrant
SMElib redesign, internal-selection batching, new line-centric opacity engine,
or `long_continuum=False` modernization was added in this release audit.

## 13. Verification and artifacts

- Full suite on the rebuilt release-source library: `187 passed, 1 xfailed,
  15 warnings` in 50.01 s.
- New release tests cover cache invalidation and adaptive transfer capacity
  sizing. Existing production tests cover immutable batched line state in both
  geometries, grazing-ray order, fixed-node equality, and spherical NLTE.
- Raw performance data: `analysis/v120_release_benchmark.json`.
- Per-run spectra/line state: `analysis/v120_release_benchmark_arrays/*.npz`.
- Deterministic accuracy summary: `analysis/v120_release_accuracy.json`.
- Audit drivers: `analysis/v120_release_benchmark.py` and
  `analysis/v120_release_accuracy.py`. Historical and production profiles are
  directly reproducible; regenerating C0B0/C1B0 requires the temporary
  controlled-baseline patch described above.
- Solar solve result: `analysis/v120_solar_wave_solve.json`.
- Supporting audits: `adaptive_continuum_grid_benchmark.md`,
  `batched_rkints_benchmark.md`, `spherical_batched_rkints_benchmark.md`,
  `rkints_second_pruning_audit.md`, and
  `adaptive_grid_range_state_audit.md`.

## 14. Release decision

1. Complete synthesis is 19.4--39.8x faster for the canonical 10 A matrix and
   17.9--43.2x faster for the Mg b/H-alpha cases.
2. Continuum caching alone is 1.5--17.2x end-to-end depending on which stage
   dominates; the canonical continuum/CDR subtask is 17.2x faster.
3. Batching alone is about 1.0x when exact continuum dominates and 2.8--3.0x
   end-to-end for line-rich cases; native line-rich transfer is 23.9x PP and
   up to 49.8x spherical/special-profile.
4. The combination has material positive interaction; naive multiplication is
   invalid.
5. Direct-synthesis continuum impact is at most `1.40e-9`; the broader CDR
   selection worst case is `7.24e-6` from a weak-line rank boundary.
6. Batching does not alter the controlled physics result: native intensities
   and final spectra are bitwise identical.
7. Historical MARK pruning, Wlim mutation, irregular-flux integration, and
   spherical ray ordering are intentional correctness changes.
8. Memory scales with line-list and transfer-node count without a cache
   explosion; the wide line-rich peak is about 4.6 GiB.
9. Cache invalidation is explicit and regression-tested.
10. The only newly exposed blocker, adaptive-grid storage capacity for a
    line-rich 800 A segment, is fixed and verified.

**GO for v1.2.0 release candidate.**  Remaining opportunities are incremental
optimization, not release blockers.

## Internal technical summary

> The performance modernization has two independent main components.
> Adaptive continuum-opacity caching reduces exact continuum evaluations by
> roughly 319--649x in the release synthesis matrix (up to about 826x in the
> CDR workload) and accelerates the canonical continuum/CDR work by 17.2x,
> with a maximum direct-synthesis normalized-flux impact of `1.40e-9` in the
> controlled matrix.  Generation-batched adaptive transfer keeps the existing
> adaptive wavelength criterion but enables sparse interval-indexed line
> evaluation, reducing candidate-line visits by 99.63--99.94% and accelerating
> line-rich native adaptive transfer by 23.91x plane-parallel and 23.88--31.35x
> spherical (39.46--49.83x in spherical Mg b/H-alpha), while remaining bitwise
> identical at controlled native transfer nodes.
