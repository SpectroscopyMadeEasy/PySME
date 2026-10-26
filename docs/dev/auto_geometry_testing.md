# Testing automatic atmosphere geometry

Automatic geometry is opt-in: set `sme.atmo.geom = None` and
`sme.atmo.depth = "RHOX"`. The structure is interpolated from the full grid;
transfer uses PP at `logg >= 3.5`, and SPH below that threshold when all
interpolation corners provide spherical geometry. See
[atmosphere configuration](../concepts/atmosphere.md#automatic-transfer-geometry).

The default remains explicit PP. This feature changes Python code only and
uses the existing SMElib interface. SPH transfer on TAU is a separate feature;
automatic mode rejects that combination.

## Fit result selection

With automatic geometry and free `logg`, bounds with `lower < 3.5 < upper` are
split into a low-gravity region ending at `nextafter(3.5, -inf)` and a
high-gravity region starting at 3.5. Other parameter bounds are retained.
Each fit starts from the supplied parameters, clipped into its region by
the existing solver. The two fits use independent model copies and Python
interpolation/wavelength caches. Their native calls run serially.

Only the successful candidate with the smallest final objective is returned,
updating the original SME object. The objective is recomputed from its final
accepted spectrum, including telluric transmission, weights, and the chosen
least-squares loss. Linear loss selects minimum chi-square; robust losses
select their corresponding objective rather than chi-square. Exactly equal
objectives retain the low-gravity candidate, which is evaluated first.
Non-finite candidates are rejected. A failed region is logged and the other
can still be returned; if both fail, solving raises an error.

This selection provides one numerical fit result. It does not establish a
unique geometry or guarantee a global minimum, and does not perform
statistical model comparison or return an alternative fit. Those analyses
belong in the calling workflow.

Intermediate candidate checkpoints are suppressed to avoid mixing the two
optimization histories. If `filename` is set, only the selected final SME
model is saved. A requested restore is applied before making the independent
candidate copies. Explicit geometry, bounds within one region, and fits
without free gravity keep the existing single-fit behavior.
Bounds ending exactly at 3.5 also retain the single-fit path: they do not
provide two regions with positive width. To compare both sides, include an
interval above 3.5 in the gravity bounds.

## Reproduce the boundary diagnostic

Use an independent checkout and Python environment, install the checkout
following the project's usual installation instructions, then run from the
repository root:

```sh
python examples/validate_auto_geometry.py --out /path/to/new/results
```

The output directory must not already exist. The default atmosphere is the
public `marcs2014.sav` grid resolved through PySME's LFS. To pin a downloaded
copy, pass `--source /path/to/grid.sav`. Run `--case 1` for only the solar
metallicity case, or `--skip-fits` for synthesis alone. `--mu` controls the
number of angular rays (default 7).

The diagnostic uses the two checked-in test line lists, with 111 lines in
total, in the 6436–6444 and 6550–6575 Å windows. These sparse line lists test
numerical behavior; they do not represent a complete observed spectrum.
It uses LTE, R = 50,000, no rotation/macroturbulence, and fixed continuum and
radial velocity. For each of three `(Teff, [M/H])` pairs it:

- Synthesizes 55 gravities around 3.5, including 3.5 ± 0.000001.
- Repeats crossings with the same synthesizer and checks that flux matches
  previous evaluations to `atol=rtol=1e-9`.
- Fits noiseless synthetic spectra at true gravities 3.45, 3.499, 3.5, 3.501,
  and 3.55, starting from both 3.3 and 3.7. Only gravity is free, with bounds
  [3.2, 3.8] and the existing 0.01 gravity Jacobian step. Uncertainties of
  1/300 supply weights; no random noise is added.
- Checks that the returned atmosphere geometry matches the accepted gravity
  and that the automatic request survives. The solver synthesizes once more
  at accepted parameters before calculating final fit results in auto mode.

Outputs include flux arrays (`*.npz`), objective curves (`curve.csv`), fit
results (`fits.csv`), geometry traces, and run metadata. The curve columns
are summed squared weighted residuals; `fits.csv` reports reduced chi-square.
A completed fit means the optimizer returned, not that it recovered the truth.

## Initial cross-boundary results

The 2026-10-26 run used the public MARCS grid with SHA256
`d8322d06b5ce30f25a53bb438578bb65299b83b2fd3f168a5b7dd80f71e1d857`
and SMElib v6.13.20. All 30 fits completed, and all repeated-crossing and
final-geometry assertions passed. Recovery was:

| Teff (K) | [M/H] | Fits with absolute gravity error > 0.003 | Maximum absolute gravity error |
| --- | --- | --- | --- |
| 4000 | +0.5 | 0 / 10 | 0.0000084 |
| 5500 | 0.0 | 2 / 10 | 0.00901 |
| 7000 | −2.5 | 4 / 10 | 0.04616 |

At 5500 K, truth 3.499 converged near 3.508 from both starts. At 7000 K,
truth 3.499 from start 3.7 converged near 3.5452; truth 3.500 and 3.501
from start 3.3 converged near 3.473. The hard transfer switch can therefore
introduce another local minimum and interfere with finite-difference steps.
These results precede the separate-region fitting described above. Spectra
retain the hard switch; each optimizer now stays in its own gravity region.
Reported derivatives and covariance retain the optimizer's Jacobian; the
final synthesis updates flux and chi-square, not the Jacobian. Uncertainties
at the boundary therefore require particular care.

## Separate-region validation

Using the same cached synthetic observations, separate-region fitting
completed all 60 candidates (two per original fit). All 30 selected results
recovered the generating geometry. Absolute gravity errors above 0.003
decreased from 6/30 to 0/30. Maximum absolute errors decreased as follows:

| Teff (K) | [M/H] | Cross-boundary maximum error | Separate-region maximum error |
| --- | --- | --- | --- |
| 4000 | +0.5 | 8.41e-6 | 6.09e-8 |
| 5500 | 0.0 | 0.009007 | 1.99e-7 |
| 7000 | −2.5 | 0.046159 | 6.20e-7 |

Atmosphere evaluations increased from 308 to 439, about 43%; this is a call
count, not a wall-clock timing measurement. The smallest difference in
summed chi-square between the two candidates was about 5.20 with the stated
weights. No exact ties occurred in this numerical sample; tie and failure
handling are covered separately by mocked regression tests.

These are consistency tests with seven rays and selected LTE models, not an
angular-convergence study or a general parameter-bias estimate. Testing with
the intended line lists, wavelength regions, NLTE setup, angular sampling,
and observed stars remains necessary before enabling auto mode in a pipeline.
