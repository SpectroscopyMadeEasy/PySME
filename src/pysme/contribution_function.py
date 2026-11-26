"""Experimental contribution functions on a prepared native transfer solution."""
from copy import deepcopy
import numpy as np


def contribution_function(sme, kind="intensity", *, quadrature=24):
    """Prepare a private SME state and return native, unbroadened contributions.

    Returns a list of segment dictionaries with wavelength (Angstrom), depth,
    depth_coordinate, mu and contributions (wavelength, ray, depth). Spherical
    results are shell integrals with depth_edges; divide by the chosen edges'
    widths to obtain a density. Line depression also returns L, Q, N and boundary.
    No observer-frame shifts, disc integration or instrumental broadening apply.
    The experimental interface may change.
    """
    if kind not in ("intensity", "line_depression"):
        raise ValueError("kind must be 'intensity' or 'line_depression'")
    from .synthesize import Synthesizer
    from .sme_synth import smelib_session
    work = deepcopy(sme)
    work.specific_intensities_only = True
    work.normalize_by_continuum = False
    work.vrad_flag = "none"
    work.cscale_flag = "none"
    synth = Synthesizer()
    results = []
    # Preparation and diagnostics share one lock and the same DLL instance.
    with smelib_session():
        for segment in range(len(work.wran)):
            synth.synthesize_spectrum(work, segments=[segment],
                                      keep_line_opacity=True,
                                      linelist_mode="all")
            dll = synth.get_dll()
            wave = np.asarray(work.wint[segment], dtype=float)
            coordinate = work.atmo.depth
            edges = np.asarray(work.atmo[coordinate], dtype=float).copy()
            spherical = str(work.atmo.geom).upper() == "SPH"
            result = dict(kind=kind, wavelength=wave.copy(), mu=np.array(work.mu),
                          depth_coordinate=coordinate,
                          intensity=np.asarray(work.sint[segment]).T.copy(),
                          continuum_intensity=np.asarray(work.cint[segment]).T.copy())
            if spherical:
                native = dll.GetLineDepressionContributionfunction(
                    wave, work.mu, len(edges), quadrature=quadrature)
                result.update(depth=(edges[:-1]+edges[1:])/2,
                              depth_edges=edges, representation="shell_integrals")
                if kind == "line_depression":
                    result.update(contributions=native.total, L=native.channels[..., 0, :],
                                  Q=native.channels[..., 4, :], N=native.channels[..., 5, :],
                                  boundary=native.boundary)
                else:
                    # Same native emission recurrence, with safe spherical folding.
                    result.update(contributions=native.channels[..., 3, :],
                                  continuum=native.channels[..., 2, :],
                                  boundary=native.metadata[..., 3],
                                  continuum_boundary=native.metadata[..., 2])
                result.update(native_metadata=native.metadata)
            elif kind == "line_depression":
                raise NotImplementedError("Line depression currently requires a spherical atmosphere")
            else:
                full, continuum = dll.GetContributionfunction(work.mu, wave)
                result.update(depth=edges, contributions=full, continuum=continuum,
                              representation="native_intensity_contribution")
            results.append(result)
    return results
