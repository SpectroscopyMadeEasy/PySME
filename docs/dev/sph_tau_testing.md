# Testing SPH+TAU transfer

This branch adds a native SMElib geometry/depth interface and Python dispatch
for spherical optical-depth transfer. It starts from the automatic-geometry
branch; default geometry remains PP. No geometry or depth default changes.

## Build and use

Use a separate checkout and environment, with Python 3.10--3.13 and the
normal [source-build prerequisites](../getting_started/installation.md#from-github-source-latest)
(C/C++ and Fortran compilers). The PySME build compiles and installs its
pinned SMElib together with the extension and native data files:

```sh
git clone --recurse-submodules --branch codex/sph-tau \
  https://github.com/SpectroscopyMadeEasy/PySME.git PySME-sph-tau
cd PySME-sph-tau
python3 -m venv .venv
. .venv/bin/activate
python -m pip install .
python -c 'from pysme.sme_synth import SME_DLL; d = SME_DLL(); print(d.SMELibraryVersion()); print(d.GetModelInputCapabilities())'
```

The submodule is pinned to SMElib commit
`e7986d63e62e348964bb3cefd0729f58d31ba580` on its `codex/sph-tau` branch.
Do not use `git submodule update --remote`: the recorded commit is the
tested pairing. Updating an existing checkout uses
`git submodule update --init --recursive` instead.

The first feature can be tested independently from PySME branch
`codex/auto-atmosphere-geometry`, in another checkout and environment. That
branch retains its old SMElib pin. This second branch contains both automatic
geometry and SPH+TAU. Neither test branch requires a develop/master merge or
a new public package release.

The capabilities must include `SPH:TAU`. The Python extension resolves the new
symbols optionally from its linked native library, so an older library can
still import and run the three historical modes. SPH+TAU against an older
library fails explicitly. Merely installing the new Python code cannot add
this capability to a downloaded old native library.

For a grid with optical depths and spherical heights:

```python
sme.atmo.geom = "SPH"  # None also works with automatic geometry
sme.atmo.interp = "TAU"
sme.atmo.depth = "TAU"
```

See [atmosphere configuration](../concepts/atmosphere.md#spherical-transfer-on-optical-depth)
for coordinate conventions and the distinction between interpolation and
transfer. This does not convert a Rosseland depth grid into a monochromatic
one. The native ABI is documented in `smelib/docs/spherical_tau.md`.

## Automated checks

From the repository root in the environment built above:

```sh
python -m pip install 'pytest>=7.2,<9'
python -m pytest -q test/test_sph_tau.py test/test_sph_tau_transfer.py \
    test/test_atmosphere_geometry.py test/test_auto_geometry_solver.py \
    test/test_atmospheres.py test/test_synthesis.py test/test_smelib_thread_safety.py
```

The new tests cover dispatch/capability checks, explicit input combinations,
old/explicit interface equivalence, mode switching, invalid reference opacity,
reference-cache stability with approximate wavelength opacity grids, raw
opacity units, continuum-scattering source consistency, unit NLTE departure
coefficients, and contribution-function sums including grazing rays.
They use the checked-in small atmosphere and two lines, without downloading
an atmosphere or NLTE grid. The existing process-global state isolation test
is expected to xfail; callers must keep complete native workflows serialized.

## Controlled coordinate consistency

```sh
python examples/validate_sph_tau.py --out /path/to/new/results
```

The output directory must not exist. The diagnostic uses the checked-in
`test/testatmo1.krz` thermodynamic profiles, smoothly resampled in log column
mass, with an imposed spherical radius of 1e11 cm and heights from 2.01e9 to
1e7 cm. This artificial geometric structure tests numerical consistency; it
is not a self-consistent MARCS giant model. Input settings are 5750 K,
logg=2.5, solar Asplund2009 abundances, microturbulence 2 km/s, and
`Ionization(0)` to preserve the supplied densities.

The 5000 Angstrom reference extinction is evaluated on an independent
288-layer grid. Its PCHIP interpolation in log mass is integrated on 30001
points to construct one shared TAU mapping. Transfer is then repeated with
24, 48, 96, 192, and 288 layers, with scattering-source mode off and on.
There are 111 fixed wavelengths from 5502.7 to 5503.8 Angstrom, two Fe/Cr
lines, and mu = [1, 0.6, 0.2, 0.05], including grazing rays.

The intensity metric is the maximum absolute RHOX/TAU difference over
wavelength divided by that ray's peak RHOX intensity, then maximized over
rays. Continuum relative errors are reported separately. Dividing a grazing
line difference by its nearly zero continuum gives a misleading measure.
The diagnostic requires the finest-grid intensity metric below 0.001,
continuum relative error below 0.002, and at least a factor 20 improvement
in intensity agreement over the coarsest grid. It saves all spectra and a
JSON report with input and native-library hashes.

## Interpretation and limits

For the 2026-10-27 local macOS arm64 build, the controlled test gave:

| Layers | Maximum intensity metric | Continuum error, source off | Continuum error, source on |
| --- | --- | --- | --- |
| 24 | 0.05697 | 0.05713 | 0.16994 |
| 48 | 0.01805 | 0.02122 | 0.05650 |
| 96 | 0.002092 | 0.002094 | 0.01007 |
| 192 | 0.0006920 | 0.0008253 | 0.003238 |
| 288 | 0.0005565 | 0.0006676 | 0.001703 |

An additional separate-process comparison used the original v6.13.20
library and the new library, the same 111 wavelengths and four rays, and
three historical modes with scattering-source mode off. Both intensity and
continuum arrays were bitwise equal in all three modes. Loading the new
extension against a copied old library reported only the three old
capabilities and rejected SPH+TAU explicitly.

These checks validate the coordinate implementation and selected compatibility
paths. They do not establish accuracy throughout MARCS parameter space,
validate line data, or require independently tabulated MARCS RHOX and TAU
spectra to match exactly. Windows/Linux builds still need their CI checks.
The inherited spherical solver returns zero for a grazing path that samples
only one radial layer (two mirrored points); this boundary behavior is
unchanged and should be considered separately when using extremely small mu.
