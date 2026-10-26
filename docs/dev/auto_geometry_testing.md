# Testing automatic atmosphere geometry

Automatic geometry is opt-in: set `sme.atmo.geom = None` and
`sme.atmo.depth = "RHOX"`. The structure is interpolated from the full grid;
transfer uses PP at `logg >= 3.5`, and SPH below that threshold when all
interpolation corners provide spherical geometry. See
[atmosphere configuration](../concepts/atmosphere.md#automatic-transfer-geometry).

The default remains explicit PP. This feature changes Python code only and
uses the existing SMElib interface. SPH transfer on TAU is a separate feature;
automatic mode rejects that combination.

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

## Initial results and limitations

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
This branch does not smooth spectra or modify the optimizer to conceal it.
Reported derivatives and covariance retain the optimizer's Jacobian; the
final synthesis updates flux and chi-square, not the Jacobian. Uncertainties
near the discontinuity therefore require particular care.

These are consistency tests with seven rays and selected LTE models, not an
angular-convergence study or a general parameter-bias estimate. Testing with
the intended line lists, wavelength regions, NLTE setup, angular sampling,
and observed stars remains necessary before enabling auto mode in a pipeline.
