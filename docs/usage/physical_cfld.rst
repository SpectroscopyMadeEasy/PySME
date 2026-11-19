Experimental physical line-depression contribution functions
==========================================================

This optional diagnostic requires the accompanying SMElib physical-CFLD branch.
Ordinary spectrum calculation does not activate it. The API is provisional.
It currently supports prepared spherical atmospheres. Plane-parallel analytical
validation uses prescribed ray fields; a general plane-parallel stellar interface
is not implemented.

Within a ``smelib_session`` encompassing ordinary model/EOS/opacity preparation,
call the same ``SME_DLL`` instance::

    result = dll.GetLineDepressionContributionfunction(wave, mu, layers)
    intensity_bins = result.total
    boundary = result.boundary
    flux_bins, flux_boundary = result.disc_flux(global_area_weights)

``wave`` uses the same wavelength convention as the native prepared state;
``mu`` contains native spherical rays, not a plane-parallel Jacobian. ``layers``
must equal the current native atmospheric node count. No model, EOS, continuum
or line-opacity preparation is performed by this method. Unsupported libraries
raise an explicit error. The method binds symbols through the already loaded
Python extension to share its linked SMElib state. Platform portability of this
binding must be reviewed; only Linux has been tested.

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
integrated depression does not establish angular/depth convergence. Current-develop stellar regression reproduces all 184 exported arrays from
eight frozen model35 jobs exactly, covering four full-feature windows and two
angular chunks each. This establishes numerical preservation for those inputs,
not independent source-physics validation or new angular/depth convergence.

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

    python -m pytest -q tests/cfld/test_prescribed.py

These tests perform no stellar synthesis. No molecular line lists, observations,
abundance fitting or velocity-analysis functionality are part of this change.
