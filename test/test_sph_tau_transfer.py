# -*- coding: utf-8 -*-
"""Native regression tests for explicit spherical TAU transfer.

The tests use the small embedded atmosphere shipped with the repository.  A
reference 5000 Angstrom optical-depth coordinate is derived from the native
continuum opacity on the same RHOX model, so RHOX and TAU describe one
controlled atmosphere rather than two independently guessed structures.
"""

from copy import deepcopy
from pathlib import Path

import numpy as np
import pytest

from pysme.abund import Abund
from pysme.atmosphere.krzfile import KrzFile
from pysme.linelist.linelist import LineList
from pysme.sme_synth import SME_DLL

from .conftest import skipif_smelib


pytestmark = [
    skipif_smelib,
    pytest.mark.filterwarnings(
        "ignore:.*EOS-computed electron density differs from the model"
    ),
]


def _make_linelist():
    lines = LineList()
    lines.add("Fe 1", 5502.9931, 0.9582, -3.047, 7.19, -6.22, 239.249)
    lines.add("Cr 2", 5503.5955, 4.1682, -2.117, 8.37, -6.49, 195.248)
    return lines


def _make_base_atmosphere():
    return KrzFile(str(Path(__file__).with_name("testatmo1.krz")))


def _make_spherical(atmo):
    result = deepcopy(atmo)
    result.geom = "SPH"
    result.radius = 1.0e11
    result.height = np.linspace(4.0e9, 0.0, len(result.rhox))
    return result


def _setup(dll, atmo, *, geom=None, depth=None, scattering=False):
    model = deepcopy(atmo)
    if geom is not None:
        model.geom = geom
    if depth is not None:
        model.depth = depth
    dll.SetLibraryPath()
    dll.InputLineList(_make_linelist())
    dll.InputModel(5750.0, 2.5 if model.geom == "SPH" else 4.5, 2.0, model)
    dll.InputAbund(Abund(monh=0.0, pattern="Asplund2009"))
    dll.Ionization(7)
    dll.SetVWscale(2.5)
    dll.SetH2broad()
    dll.InputWaveRange(5502.5, 5504.2)
    dll.SetContinuumScatteringSourceMode(scattering)
    dll.Opacity()
    return model


def _spectrum(dll):
    nw, wave, intensity, continuum = dll.Transf(
        [1.0], accrt=0.1, accwi=0.1
    )
    return {
        "nw": nw,
        "wave": np.asarray(wave).copy(),
        "intensity": np.asarray(intensity).copy(),
        "continuum": np.asarray(continuum).copy(),
    }


def _assert_spectrum_equal(actual, expected):
    assert actual["nw"] == expected["nw"]
    np.testing.assert_array_equal(actual["wave"], expected["wave"])
    np.testing.assert_allclose(
        actual["intensity"], expected["intensity"], rtol=2e-12, atol=2e-14
    )
    np.testing.assert_allclose(
        actual["continuum"], expected["continuum"], rtol=2e-12, atol=2e-14
    )


@pytest.fixture(scope="module")
def tau_atmosphere():
    """Return a base PP/RHOX model and its native 5000 A TAU coordinate."""

    base = _make_base_atmosphere()
    dll = SME_DLL()
    model = _setup(dll, base, geom="PP", depth="RHOX")
    _, _, chi = dll.GetContinuumOpacityComponents(5000.0)
    chi = np.asarray(chi, dtype=float)
    rhox = np.asarray(model.rhox, dtype=float)
    assert np.all(np.isfinite(chi))
    assert np.all(chi > 0)
    assert np.all(np.diff(rhox) > 0)

    # Integral of chi dm, with a small positive offset so the first point is
    # strictly inside the native grid instead of exactly tau=0.
    tau = np.empty_like(rhox)
    tau[0] = chi[0] * rhox[0] * 1e-6
    tau[1:] = tau[0] + np.cumsum(
        0.5 * (chi[1:] + chi[:-1]) * np.diff(rhox)
    )
    assert np.all(np.isfinite(tau))
    assert np.all(tau > 0)
    assert np.all(np.diff(tau) > 0)

    tau_model = deepcopy(base)
    tau_model.tau = tau
    tau_model.depth = "TAU"
    return base, tau_model, chi


@pytest.mark.parametrize(
    "geom,depth",
    [("PP", "RHOX"), ("PP", "TAU"), ("SPH", "RHOX")],
)
def test_legacy_and_explicit_input_match_for_old_modes(
    tau_atmosphere, geom, depth
):
    base, tau_model, _ = tau_atmosphere
    source = tau_model if depth == "TAU" else base
    if geom == "SPH":
        source = _make_spherical(source)

    legacy = SME_DLL()
    _setup(legacy, source, geom=geom, depth=depth)
    expected = _spectrum(legacy)

    explicit = SME_DLL()
    model = deepcopy(source)
    model.geom = geom
    model.depth = depth
    # InputModel dispatches SPH+TAU only; for the three historical modes call
    # the explicit method directly to compare the native Ex ABI.
    explicit.SetLibraryPath()
    explicit.InputLineList(_make_linelist())
    explicit.InputModelEx(
        5750.0,
        2.5 if geom == "SPH" else 4.5,
        2.0,
        model,
    )
    explicit.InputAbund(Abund(monh=0.0, pattern="Asplund2009"))
    explicit.Ionization(7)
    explicit.SetVWscale(2.5)
    explicit.SetH2broad()
    explicit.InputWaveRange(5502.5, 5504.2)
    explicit.SetContinuumScatteringSourceMode(False)
    explicit.Opacity()
    actual = _spectrum(explicit)
    _assert_spectrum_equal(actual, expected)


def test_sph_tau_explicit_transfer_is_finite(tau_atmosphere):
    _, tau_model, _ = tau_atmosphere
    sph_tau = _make_spherical(tau_model)
    dll = SME_DLL()
    _setup(dll, sph_tau, geom="SPH", depth="TAU")
    result = _spectrum(dll)

    assert np.all(np.isfinite(result["wave"]))
    assert np.all(np.isfinite(result["intensity"]))
    assert np.all(np.isfinite(result["continuum"]))
    assert np.all(result["continuum"] > 0)


def test_sph_tau_rhox_mode_switch_is_stable(tau_atmosphere):
    base, tau_model, _ = tau_atmosphere
    dll = SME_DLL()
    sequence = (
        ("SPH", "TAU", _make_spherical(tau_model)),
        ("PP", "RHOX", base),
        ("SPH", "RHOX", _make_spherical(base)),
        ("SPH", "TAU", _make_spherical(tau_model)),
    )
    observations = []
    for geom, depth, model in sequence:
        _setup(dll, model, geom=geom, depth=depth)
        observations.append(_spectrum(dll))

    _assert_spectrum_equal(observations[0], observations[-1])
    repeat = SME_DLL()
    _setup(repeat, _make_spherical(tau_model), geom="SPH", depth="TAU")
    _assert_spectrum_equal(observations[-1], _spectrum(repeat))


@pytest.mark.parametrize("depth", ["RHOX", "TAU"])
@pytest.mark.parametrize("mu", [0.2, 0.05])
def test_sph_grazing_rays_are_finite(tau_atmosphere, depth, mu):
    base, tau_model, _ = tau_atmosphere
    source = _make_spherical(tau_model if depth == "TAU" else base)
    dll = SME_DLL()
    _setup(dll, source, geom="SPH", depth=depth)
    nw, wave, intensity, continuum = dll.Transf(
        [1.0, mu], accrt=0.1, accwi=0.1
    )

    assert nw == len(wave)
    assert np.all(np.isfinite(wave))
    assert np.all(np.isfinite(intensity))
    assert np.all(np.isfinite(continuum))


def test_sph_tau_rejects_missing_reference_opacity():
    model = _make_spherical(_make_base_atmosphere())
    model.depth = "TAU"
    model.tau = np.geomspace(1e-5, 1.0, len(model.rhox))
    model.opflag = np.zeros(20, dtype=np.int16)
    dll = SME_DLL()
    dll.SetLibraryPath()
    dll.InputLineList(_make_linelist())
    dll.InputModel(5750.0, 2.5, 2.0, model)
    dll.InputAbund(Abund(monh=0.0, pattern="Asplund2009"))
    dll.Ionization(7)
    dll.InputWaveRange(5502.5, 5504.2)

    with pytest.raises((RuntimeError, ValueError), match="[Ii]nvalid.*reference"):
        dll.Opacity()
    with pytest.raises((RuntimeError, ValueError), match="[Ii]nvalid.*reference"):
        dll.GetContinuumScatteringSource(5000.0)


def test_sph_tau_reference_opacity_cache_survives_transfer(tau_atmosphere):
    _, tau_model, _ = tau_atmosphere
    dll = SME_DLL()
    _setup(dll, _make_spherical(tau_model), geom="SPH", depth="TAU")
    if not hasattr(dll, "SetContinuumOpacityGrid"):
        pytest.skip("continuum opacity grid controls are unavailable")
    dll.SetContinuumOpacityGrid("adaptive", rtol=1e-3)
    before = np.asarray(dll.GetOpacity("COPSTD")).copy()
    _spectrum(dll)
    after = np.asarray(dll.GetOpacity("COPSTD")).copy()
    np.testing.assert_allclose(after, before, rtol=1e-12, atol=0.0)


def test_sph_tau_line_opacity_matches_rhox_for_same_structure(tau_atmosphere):
    base, tau_model, _ = tau_atmosphere
    rhox = SME_DLL()
    _setup(rhox, _make_spherical(base), geom="SPH", depth="RHOX")
    rhox.GetLineOpacity(5502.9931)
    rhox_lop, rhox_cop, rhox_scr, rhox_tsf, rhox_csf = rhox.GetLineOpacity(
        5502.9931
    )

    tau = SME_DLL()
    _setup(tau, _make_spherical(tau_model), geom="SPH", depth="TAU")
    tau_lop, tau_cop, tau_scr, tau_tsf, tau_csf = tau.GetLineOpacity(
        5502.9931
    )

    for actual, expected in (
        (tau_lop, rhox_lop),
        (tau_cop, rhox_cop),
        (tau_scr, rhox_scr),
        (tau_tsf, rhox_tsf),
        (tau_csf, rhox_csf),
    ):
        np.testing.assert_allclose(actual, expected, rtol=3e-3, atol=1e-12)


def test_sph_tau_scattering_source_matches_line_opacity_units(tau_atmosphere):
    _, tau_model, _ = tau_atmosphere
    dll = SME_DLL()
    _setup(
        dll,
        _make_spherical(tau_model),
        geom="SPH",
        depth="TAU",
        scattering=True,
    )
    line_wave = 5502.9931
    dll.Transf([1.0], wave=np.array([line_wave]), accrt=0.1, accwi=0.1)
    jbar, source = dll.GetContinuumScatteringSource(line_wave)
    lop, cop, scr, tsf, csf = dll.GetLineOpacity(line_wave)
    kappa, sigma, chi = dll.GetContinuumOpacityComponents(line_wave)

    np.testing.assert_allclose(csf, source, rtol=2e-12, atol=2e-14)
    np.testing.assert_allclose(cop, chi, rtol=2e-12, atol=2e-14)
    np.testing.assert_allclose(scr, sigma, rtol=2e-12, atol=2e-14)
    assert np.all(np.isfinite(jbar))
    assert np.all(np.isfinite(lop))
    assert np.all(cop > 0)
    # These are raw opacity/source values, not continuum-normalized profiles.
    assert np.max(cop) > 1.0
    assert np.max(source) > 1.0


def test_sph_tau_contribution_function_integrates_to_specific_intensity(
    tau_atmosphere,
):
    _, tau_model, _ = tau_atmosphere
    dll = SME_DLL()
    _setup(dll, _make_spherical(tau_model), geom="SPH", depth="TAU")
    mu = np.array([1.0, 0.2, 0.05])
    wave = np.array([5502.9931, 5503.5955])
    _, _, intensity, continuum = dll.Transf(
        mu, wave=wave, accrt=0.1, accwi=0.1
    )

    contribution, continuum_contribution = dll.GetContributionfunction(
        mu, wave, long_continuum=True
    )
    assert contribution.shape == (wave.size, mu.size, len(tau_model.rhox))
    assert continuum_contribution.shape == contribution.shape
    assert np.all(np.isfinite(contribution))
    assert np.all(np.isfinite(continuum_contribution))
    np.testing.assert_allclose(
        contribution.sum(axis=2), intensity.T, rtol=2e-12, atol=1e-10
    )
    np.testing.assert_allclose(
        continuum_contribution.sum(axis=2),
        continuum.T,
        rtol=2e-12,
        atol=1e-10,
    )


def test_sph_tau_unit_departures_preserve_lte_transfer(tau_atmosphere):
    _, tau_model, _ = tau_atmosphere
    dll = SME_DLL()
    try:
        _setup(dll, _make_spherical(tau_model), geom="SPH", depth="TAU")
        wave = np.array([5502.9931, 5503.5955])
        _, _, lte_intensity, lte_continuum = dll.Transf(
            [1.0], wave=wave, accrt=0.1, accwi=0.1
        )
        dll.InputDepartureCoefficients(np.ones((len(tau_model.rhox), 2)), 0)
        _, _, nlte_intensity, nlte_continuum = dll.Transf(
            [1.0], wave=wave, accrt=0.1, accwi=0.1
        )
        np.testing.assert_allclose(
            nlte_intensity, lte_intensity, rtol=3e-12, atol=2e-10
        )
        np.testing.assert_allclose(
            nlte_continuum, lte_continuum, rtol=3e-12, atol=2e-10
        )
    finally:
        dll.ResetDepartureCoefficients()
        dll.SetContinuumScatteringSourceMode(False)


@pytest.mark.parametrize("scattering", [False, True])
def test_sph_tau_transfer_is_finite_with_scattering(tau_atmosphere, scattering):
    _, tau_model, _ = tau_atmosphere
    dll = SME_DLL()
    _setup(
        dll,
        _make_spherical(tau_model),
        geom="SPH",
        depth="TAU",
        scattering=scattering,
    )
    result = _spectrum(dll)
    assert np.all(np.isfinite(result["intensity"]))
    assert np.all(np.isfinite(result["continuum"]))
