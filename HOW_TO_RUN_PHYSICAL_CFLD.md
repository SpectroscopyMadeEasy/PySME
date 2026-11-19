# How to run the experimental physical CFLD

This guide describes the review branch, not a released PySME feature. You need **both the modified PySME Python code and its companion modified SMElib**. Installing ordinary PySME alone does not supply this API.

## Existing emission CF versus the new line-depression CF

These APIs answer different scientific questions. The new method supplements, rather than replaces, `GetContributionfunction`.

| Aspect | Existing `GetContributionfunction(mu, wave, ...)` | New `GetLineDepressionContributionfunction(wave, mu, layers, ...)` |
|---|---|---|
|Question|Which atmospheric regions contribute to the emergent radiation?|Which regions contribute to the absolute continuum-minus-full absorption depression?|
|Quantity|Emission contribution decomposition for full-opacity intensity and, when requested, continuum intensity|Transfer-equation decomposition of Ic−If, with disc integration giving Fc−Ff|
|Typical purpose|Study the depth weighting of escaping radiation at a wavelength|Study the depth weighting of observable absorption relative to the model continuum|
|Return contract|Full table and continuum table when `long_continuum=True`; full table alone otherwise|Structured result: six channels, complete L+Q+N shell integrals and separate boundary/metadata|
|Key local information|Native emission/source contributions and attenuation|Internal directional continuum intensity, line extinction, native source differences and attenuation, with numerical interpolation accounting|
|Status|Existing PySME API, unchanged by this patch|Experimental API requiring the accompanying modified native library|

The call signatures have different argument order: **existing: mu, wave; new: wave, mu, layers**. The high-level synthesis `contribution_function` option still requests the existing emission diagnostic; it does not activate the new physical CFLD.

Subtracting compatible absolute continuum and full emission contributions can close the total depression **when the corresponding boundaries are included**, but does not establish the same depth-local distribution as the physical line-depression equation. Its boundary difference also differs from the physical depression boundary, which uses full-opacity attenuation. Do not directly subtract returned legacy arrays without auditing their normalisation, coordinate measure, boundary handling and geometry; the existing array conventions are not automatically those of the new absolute shell integrals.

The new result's emission-control channels are extracted from the traced native recurrence. They are not obtained by calling the old emission-CF interface. This avoids reliance on that interface's previously identified spherical path-length/output-stride issue; it does not claim to repair or independently revalidate the old interface.

“Physical” identifies the continuum-minus-full driving definition. The implemented total includes numerical correction N; N is not an additional physical process or an uncertainty estimate. Closure and exact port regression do not establish full source-physics consistency or depth convergence. Neither diagnostic alone measures the spatial abundance of a species, identifies distinct atmospheric shells, or determines gas velocities. Peak depths from the two definitions need not coincide.

## 1. Choose the appropriate starting point

| Task | Required inputs | Starts stellar synthesis? |
|---|---|---|
|Prescribed-field tests|Compiled diagnostic library, NumPy, pytest|No|
|Analytical benchmark|Same library; benchmark harness and SciPy|No stellar synthesis; prescribed transfer only|
|Stellar CFLD|Prepared spherical atmosphere/EOS/opacities and native rays in the same library state|Yes, preparation is required; the diagnostic repeats native formal transfer|
|Review saved regression|Existing outputs and local audit script|No|

The prototype has been built and tested on Linux. Native export/build integration and cross-platform symbol lookup remain under maintainer review. The local build scripts are audit recipes, not a portable installation procedure. Do not replace your working installed library to try the diagnostic.

## 2. Run the prescribed-field tests first

From this PySME checkout, with an existing isolated diagnostic library:

```bash
CFLD_LIBRARY=/absolute/path/to/libsme_cfld.so python -m pytest -q tests/cfld/test_prescribed.py
```

Expected result for the reviewed version: **22 passed**. Without CFLD_LIBRARY, native tests skip; a test run containing skips is not equivalent to this validation. This command does not initialise a stellar atmosphere or call EOS/opacity preparation.

The tests load the local Python interface directly. For normal imports in your own scripts, ensure this branch's `src` directory is selected ahead of the installed PySME. Confirm `pysme.__file__` and the loaded native bridge identity; a displayed library version alone is insufficient evidence that the intended bridge/state is in use.

## 3. Minimal prescribed-path example

With this branch importable and the isolated library available:

```python
import numpy as np
from pysme.line_depression import NativeLineDepression

native = NativeLineDepression('/absolute/path/to/libsme_cfld.so')
m = np.geomspace(0.001, 20.0, 56)
kc = np.ones_like(m)
kf = 2.0 * kc
B = 1.0 + 0.5 * m
result = native.prescribed_path(m, kc, kf, B, B, B, quadrature=24)

print('Native depression:', result.metadata[0] - result.metadata[1])
print('Integrated diagnostic:', result.total.sum() + result.boundary)
```

These are controlled mathematical fields, **not a stellar atmosphere or synthetic observation**. Extinction times path-column-mass increments gives optical thickness; native output scaling is retained. The inherited outgoing interval 0–1 omission also applies. The example checks the mechanics; the independent analytical benchmark tests a specified physical limit.

## 4. Request a stellar CFLD from a prepared native state

The following is an **API usage fragment**, not a complete stellar setup. `dll`, `sme`, `wave`, `mu`, `global_area_weights` and `logtauRoss` must come from a consistent caller-owned setup. The compiled Python bridge must link to the diagnostic SMElib; a second unrelated library state will not work.

```python
import numpy as np
from pysme.sme_synth import smelib_session
from pysme.synthesize import Synthesizer

with smelib_session():
    # sme already contains the atmosphere, absolute abundances, line list,
    # wavelength grid, ray directions and explicit synthesis settings.
    # dll is the verified diagnostic-enabled instance.
    synthesizer = Synthesizer(dll=dll)
    returned = synthesizer.synthesize_spectrum(
        sme, linelist_mode='all', contribution_function=False
    )
    if hasattr(returned, 'wave'):
        sme = returned

    ordinary = dll.Transf(
        mu, wave=wave, keep_lineop=True, long_continuum=True
    )
    result = dll.GetLineDepressionContributionfunction(
        wave, mu, layers=len(logtauRoss), quadrature=24
    )

flux_bins, flux_boundary = result.disc_flux(global_area_weights)
cfld_per_dex_ross = flux_bins / np.abs(np.diff(logtauRoss))
continuum_flux = np.pi * np.einsum('r,rw->w', global_area_weights,
                                  ordinary[3])
full_flux = np.pi * np.einsum('r,rw->w', global_area_weights, ordinary[2])
closure_residual = flux_bins.sum(axis=-1) + flux_boundary - (
    continuum_flux - full_flux
)
```

Use the same wavelength grid and native rays in all calls. These are **absolute, unbroadened** quantities; the local validated runner explicitly disables continuum normalisation, RV fitting and broadening. Do not insert a normalised or convolved spectrum into the closure comparison. This fragment does not certify arbitrary setup compatibility.

The diagnostic requires a prepared spherical state and correct node count. It evaluates native opacity/source fields and repeats the formal solution with tracing enabled; it is not a function of the emergent spectrum alone. Native status 10 before preparation is an expected rejection, not a physical-model result. Execute complete workflows under the process lock; thread-local tracing does not protect global atmospheric state.

## 5. Understand and save the outputs

- `result.channels`: `(nwave, nray, 6, nradial-1)`; channels L, Q+N, continuum emission, full emission, Q, N.
- `result.total`: signed L+Q+N shell integrals per ray, not a density.
- `result.boundary`: separate physical depression boundary per ray.
- `result.metadata`: native intensity/boundary/internal-minimum/path metadata; see the reviewer guide.
- `result.disc_flux(weights)`: pi-weighted flux-bin and boundary arrays.

If rays must be chunked (the tested native maximum is 77 per call), retain **global** area weights and sum chunk fluxes. Never normalise each chunk to one. For log tau5000 density, divide by its own widths rather than relabelling the Rosseland density. Keep boundaries separate and preserve negative contributions.

Save wavelength convention, arrays, coordinates, ray nodes/weights, source mode, atmosphere/abundance identities, native/bridge/Python versions and hashes, and ordinary spectra. Verify ray and disc closure and native ordinary-before/after agreement. Assess distribution convergence separately. Do not rescale CFLDs to close the integral.

## 6. Limitations to retain in any result

The outer interval omission and two-entry grazing behaviour are preserved. Full and continuum source treatment need not be equivalent self-consistent scattering solutions; general Q is retained. N is local numerical interpolation accounting, not an uncertainty. The current exact stellar regression demonstrates preservation for frozen inputs, not universal physical validation or depth convergence. Native function names and public packaging remain provisional.

See `docs/usage/physical_cfld.rst` and the companion SMElib `PHYSICAL_CFLD_REVIEWER_GUIDE.md` for the method and review questions.
