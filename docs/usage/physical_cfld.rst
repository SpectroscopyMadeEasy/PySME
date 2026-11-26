Experimental physical line-depression contribution functions
==========================================================

This optional diagnostic requires the accompanying SMElib physical-CFLD branch.
Ordinary spectrum calculation does not activate it. The API is provisional.
It currently supports prepared spherical atmospheres. Plane-parallel analytical
validation uses prescribed ray fields; a general plane-parallel stellar interface
is not implemented.

The implementation aims to follow the absolute line-depression contribution
approach of `Kravchenko et al. (2018), A&A 610, A29
<https://doi.org/10.1051/0004-6361/201731530>`_. It independently adapts that
physical objective to SMElib's actual source functions, spherical transfer and
numerical interpolation. It does not claim identical source assumptions or
complete physical validation; the native numerical correction N is distinct
from the underlying physical formulation.

Minimal executable example
--------------------------

After installing the modified Python package and building its companion native
library, run from the PySME checkout::

    python examples/line_depression.py --library /absolute/path/to/libsme_cfld.so

Alternatively set ``CFLD_LIBRARY``. Add ``--output demo_cfld.npz`` to save the
six channels, boundary metadata, complete interval contributions and reference-
optical-depth density. Existing output files are refused. Expected output is
``PASS_PRESCRIBED_DEMO``, a relative closure residual below 1e-6, and a passing
no-line cancellation check. The tested Linux build gives a residual about
-4.53e-9; last digits need not be identical on other builds.

This deliberately simple single-ray linear-source example needs no atmosphere,
line list, downloaded stellar data or EOS calculation. Its depth coordinate is
reference continuum optical depth, not tauRoss or tau5000. It demonstrates the
native diagnostic and accounting, not a stellar spectrum or full analytical
shape-validation suite. A real stellar setup requires the preparation below.

Prerequisites and calling the diagnostic
----------------------------------------

Use the modified PySME code and a Python bridge linked to the companion
diagnostic-enabled SMElib. An unmodified installed library does not provide
these symbols. Do not replace a working installation to test this prototype.

The following is an API fragment, not a complete stellar setup. ``dll`` must
be the same verified instance used to prepare the atmosphere, abundances,
line list, EOS and opacities. Keep continuum normalisation, broadening and
RV fitting disabled when testing absolute native closure.

Within a ``smelib_session`` encompassing ordinary model/EOS/opacity preparation,
call the same ``SME_DLL`` instance::

    result = dll.GetLineDepressionContributionfunction(wave, mu, layers)
    intensity_bins = result.total
    boundary = result.boundary
    flux_bins, flux_boundary = result.disc_flux(global_area_weights)
    cfld_ross = flux_bins / abs(np.diff(logtauRoss))

Here ``np`` is NumPy. The required closure is ``flux_bins.sum(axis=-1) +
flux_boundary == continuum_flux - full_flux`` within the declared numerical
tolerance, using matching absolute, unbroadened native fluxes.

``wave`` uses the same wavelength convention as the native prepared state;
``mu`` contains native spherical rays, not a plane-parallel Jacobian. ``layers``
must equal the current native atmospheric node count. No model, EOS, continuum
or line-opacity preparation is performed by this method. Unsupported libraries
raise an explicit error. The method binds symbols through the already loaded
Python extension to share its linked SMElib state. Platform portability of this
binding must be reviewed; only Linux has been tested.

The diagnostic evaluates opacity/source fields and repeats native formal
transfer with tracing enabled; it is not an accessor for previously cached
intensity traces. A request before model/EOS/opacity preparation fails explicitly.

Outputs and conventions
-----------------------

Raw channels have shape (wavelength, ray, 6, radial_interval). The channels are
line driving L, Q+N, continuum emission, full emission, Q, N. Complete physical
contributions are L+Q+N; emission-continuum minus emission-full is a different
local decomposition. Boundary contributions are never assigned to a shell.

The general source term is q = chi_c S_c - chi_f S_f + chi_l B. It is not assumed
to equal a difference of scattering mean intensities. L uses the internal
continuum intensity and thermal-source contrast. Q is the nodal general-source
term; N accounts locally for the native Bezier interpolation relative to the
chosen nodal representation. N is neither a fitted closure residual nor an
uncertainty estimate. Signed values are retained.

Native output retains the solver's FLUX_SCALE conversion and absolute specific
intensity convention. The result is not continuum normalised. Disc integration
uses pi times global projected-area weights and reports flux at the outer model
radius; caller-supplied chunks must retain their global weights. Coordinate
conversion divides shell integrals by each coordinate's own positive widths.
For flux density per dex, divide flux_bins by abs(diff(logtauRoss)) or separately
abs(diff(logtau5000)); boundaries remain separate. No broadening is applied.

The existing native transfer domain excludes interval 0--1 and assigns zero
output to two-entry grazing paths. This implementation preserves those
conventions. Near/far paths are folded onto radial intervals. Agreement of
integrated depression does not establish angular/depth convergence.
Compared with the existing emission contribution function
--------------------------------------------------------

``GetContributionfunction(mu, wave, ...)`` is the existing emission diagnostic:
it describes contributions to emergent full-opacity radiation and, with
``long_continuum=True``, continuum radiation. The new
``GetLineDepressionContributionfunction(wave, mu, layers, ...)`` describes
contributions to absolute continuum-minus-full line depression. Notice the
different argument order. The high-level ``contribution_function`` option still
selects the existing emission diagnostic; the new API must be requested explicitly.

Use the emission CF to investigate escaping radiation, and the line-depression
CF to investigate absorption weighting relative to the model continuum. An
absolute difference of emission decompositions may recover the total depression
with compatible normalisation and boundary terms, but its local contributions
need not equal the physical depression-driving decomposition. Do not subtract
legacy returned tables without checking their own units, normalisation,
coordinates and boundaries.

The new emission-control channels come from the traced recurrence, not the old
public CF routine. The old API is unchanged, including its outstanding spherical
indexing concern. Neither CF directly gives a species abundance distribution or
a gas velocity. “Physical” specifies the driving definition; the complete native
result also retains the numerical interpolation correction N. It is not a claim
of universal physical or numerical convergence validation.

Portable prescribed-opacity tests
---------------------------------

Set CFLD_LIBRARY to an isolated compiled diagnostic library and run::

    python -m pytest -q test/test_line_depression.py

These tests perform no stellar synthesis. No molecular line lists, observations,
abundance fitting or velocity-analysis functionality are part of this change.

Save wavelength convention, atmosphere/abundance identities, ray nodes/weights,
coordinate arrays, source mode, native/bridge versions and ordinary spectra
alongside the diagnostic. Preserve signed values and do not rescale to force
closure. Native ray limits require chunking when exceeded; the tested limit is
77 rays per call. Keep global weights when summing chunks.

The companion SMElib document ``docs/physical_cfld.md`` contains the derivation,
L/Q/N definitions, spherical mapping and boundary conventions. Detailed
validation history is supplied with the PR rather than maintained in this page.

The test lives in the existing ``test/`` directory and is discovered by the
repository's pytest configuration. The standard command loads the shared
``test/conftest.py`` and therefore requires PySME's normal test dependencies.
For a minimal prescribed-field check that bypasses those unrelated fixtures,
use ``python -m pytest --noconftest -q test/test_line_depression.py`` with
``CFLD_LIBRARY`` set as above. This restricted run does not validate integration
with the shared test fixtures.

Experimental common entry point
-------------------------------

This interface is experimental and may change before release::

    from pysme.contribution_function import contribution_function
    segments = contribution_function(sme, kind="line_depression")
    emission = contribution_function(sme, kind="intensity")

The function prepares a private copy of the SME structure under one native
session lock. Each segment returns wavelength, depth, mu and contributions.
Spherical arrays have shape (wavelength, ray, radial interval); depth_edges
identifies their integration bins. They are shell integrals, not per-dex
densities. Line depression additionally returns L, Q, N and boundary terms.
Native wavelength and intensity units are retained, before broadening or RV
alignment. The input SME structure is not modified. Spherical intensity uses
the same native emission recurrence with safe radial-shell folding; the legacy
GetContributionfunction method is unchanged. Plane-parallel intensity retains
its native depth-node representation; plane-parallel line depression is not
supported by the atmospheric diagnostic.
