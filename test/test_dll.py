# -*- coding: utf-8 -*-
from copy import deepcopy
from pathlib import Path
from os.path import dirname

import numpy as np
import pytest

from pysme.abund import Abund
from pysme.atmosphere.krzfile import KrzFile, atmoic_mass
from pysme.linelist.linelist import LineList
from pysme.sme_synth import SME_DLL


# Create Objects to pass to library
# Their functionality is tested in other test files, so we assume it works
@pytest.fixture
def cwd():
    return dirname(__file__)


@pytest.fixture
def libsme():
    dll = SME_DLL()
    # SMElib keeps this mode in process-global state, so isolate tests that
    # intentionally enable strict precomputed-line-info handling.
    dll.SetLineInfoMode(0)
    return dll


@pytest.fixture
def teff():
    return 5750


@pytest.fixture
def grav():
    return 4.5


@pytest.fixture
def vturb():
    return 2


@pytest.fixture
def wfirst():
    return 5500


@pytest.fixture
def wlast():
    return 5600


@pytest.fixture
def vw_scale():
    return 2.5


@pytest.fixture
def mu():
    return [1]


@pytest.fixture
def accrt():
    return 0.1


@pytest.fixture
def accwt():
    return 0.1


@pytest.fixture
def linelist():
    linelist = LineList()
    linelist.add("Fe 1", 5502.9931, 0.9582, -3.047, 7.19, -6.22, 239.249)
    linelist.add("Cr 2", 5503.5955, 4.1682, -2.117, 8.37, -6.49, 195.248)
    linelist.add("Mn 1", 5504.0000, 2.0000, -3.000, 8.00, -6.50, 200.247)
    return linelist


@pytest.fixture
def atmo(cwd):
    return KrzFile(cwd + "/testatmo1.krz")


@pytest.fixture
def abund():
    return Abund(monh=0, pattern="Asplund2009")


def test_basic(libsme, wfirst, wlast, vw_scale):
    """Test instantiation of library object and some basic functions"""
    print(libsme.SMELibraryVersion())
    libsme.InputWaveRange(wfirst, wlast)
    libsme.SetVWscale(vw_scale)
    libsme.SetH2broad()

    # assert libsme.file
    # assert libsme.wfirst == wfirst
    # assert libsme.wlast == wlast
    # assert libsme.vw_scale == vw_scale
    # assert libsme.h2broad


def test_eos_warm_start_mode_api(libsme):
    assert libsme.SetEosWarmStartMode(True)
    assert libsme.SetEosWarmStartMode(False)


def test_eos_warm_start_mode_tolerates_older_extension(libsme, monkeypatch):
    import pysme.sme_synth as sme_synth

    monkeypatch.setattr(sme_synth._smelib, "SetEosWarmStartMode", None)
    assert not libsme.SetEosWarmStartMode(True)


def test_select_strong_lines_by_bins_native_api():
    wavelength = np.array([5000.0, 5000.01, 5000.02, 5000.2, 5000.21])
    metric = np.array([2e-4, 3e-4, 8e-4, 4e-4, np.nan])
    valid = np.array([True, True, True, True, True])

    strong = SME_DLL.SelectStrongLinesByBins(
        wavelength,
        metric,
        bin_width=0.2,
        threshold=5e-4,
        valid_mask=valid,
    )

    assert np.array_equal(strong, [False, False, True, False, False])
    assert SME_DLL.SelectStrongLinesByBins([], []).size == 0
    with pytest.raises(RuntimeError, match="bin_width"):
        SME_DLL.SelectStrongLinesByBins(wavelength, metric, bin_width=0)


def test_linelist(libsme, linelist):
    """Test linelist behaviour"""
    libsme.InputLineList(linelist)
    outlist = libsme.OutputLineList()

    for i in range(len(linelist)):
        outline = [x for x in outlist[i]]
        outline[1] = np.log10(outline[1])

    # TODO
    # libsme.UpdateLineList()

    with pytest.raises(TypeError):
        libsme.InputLineList(None)


def test_atmosphere(libsme, atmo, teff, grav, vturb):
    """Test atmosphere behaviour"""

    libsme.InputModel(teff, grav, vturb, atmo)

    # TODO test different geometries
    with pytest.raises(ValueError):
        libsme.InputModel(-1000, grav, vturb, atmo)

    with pytest.raises(ValueError):
        libsme.InputModel(teff, grav, -10, atmo)

    with pytest.raises(TypeError):
        libsme.InputModel(teff, grav, vturb, None)


def test_abund(libsme, abund):
    """Test abundance behaviour"""
    libsme.InputAbund(abund)

    # TODO: What should be the expected behaviour?
    empty = Abund(monh=0, pattern="empty")
    empty.update_pattern({"H": 12})
    libsme.InputAbund(empty)

    with pytest.raises(TypeError):
        libsme.InputAbund(None)


def test_krzfile_mu_uses_number_ratios(atmo):
    ratios = atmo.abund.get_pattern(type="n/nH")
    valid = [(element, value) for element, value in ratios.items() if not np.isnan(value)]
    expected_mu = sum(
        value * atmoic_mass[element] for element, value in valid
    ) / sum(value for _, value in valid)

    assert np.isclose(atmo.get_mu_from_abund(), expected_mu)


def test_krzfile_reads_monh_from_abundance_scale(tmp_path, cwd, atmo):
    source = Path(cwd) / "testatmo1.krz"
    content = source.read_text()
    content = content.replace("TITLE  [0.0]", "TITLE       ", 1)
    content = content.replace("ABUNDANCE SCALE   1.00000", "ABUNDANCE SCALE   0.10000", 1)
    target = tmp_path / "scaled_from_abundance.krz"
    target.write_text(content)

    scaled = KrzFile(str(target))

    assert np.isclose(scaled.monh, -1.0)
    assert np.isclose(scaled.abundance_scale, 0.1)
    assert np.isclose(scaled.abund.A["Fe"], atmo.abund.A["Fe"] - 1.0)


def test_krzfile_warns_when_header_and_abundance_scale_disagree(tmp_path, cwd):
    source = Path(cwd) / "testatmo1.krz"
    content = source.read_text()
    content = content.replace("ABUNDANCE SCALE   1.00000", "ABUNDANCE SCALE   0.10000", 1)
    target = tmp_path / "mismatch_scale_header.krz"
    target.write_text(content)

    with pytest.warns(UserWarning, match="ATLAS abundance scale and header metallicity disagree"):
        atmo = KrzFile(str(target))

    assert np.isclose(atmo.monh, -1.0)


def test_transf(
    libsme,
    linelist,
    teff,
    grav,
    vturb,
    atmo,
    abund,
    vw_scale,
    wfirst,
    wlast,
    mu,
    accrt,
    accwt,
):
    """Test radiative transfer"""
    libsme.SetLibraryPath()

    libsme.InputLineList(linelist)
    libsme.InputModel(teff, grav, vturb, atmo)
    libsme.InputAbund(abund)
    libsme.Ionization(0)
    libsme.SetVWscale(vw_scale)
    libsme.SetH2broad()

    libsme.InputWaveRange(wfirst, wlast)
    libsme.Opacity()

    nw, wave, synth, cont = libsme.Transf(mu, accrt=accrt, accwi=accwt)
    assert nw == len(wave) == synth.shape[-1] == cont.shape[-1]
    assert nw > 0
    assert np.isclose(wave[0], wfirst)
    assert np.isclose(wave[-1], wlast)
    assert np.all(np.diff(wave) > 0)

    density = libsme.GetDensity()
    assert np.allclose(density, atmo.rho, rtol=3e-1, equal_nan=True)

    xne = libsme.GetNelec()
    assert np.allclose(xne, atmo.xne, rtol=2e-1)

    xna = libsme.GetNatom()
    assert np.allclose(xna, atmo.xna, rtol=1e-1)

    lop, cop, scr, tsf, csf = libsme.GetLineOpacity(linelist.wlcent[0])
    kappa, sigma, chi = libsme.GetContinuumOpacityComponents(linelist.wlcent[0])
    assert kappa.shape == sigma.shape == chi.shape == cop.shape == scr.shape
    assert np.all(np.isfinite(kappa))
    assert np.all(np.isfinite(sigma))
    assert np.all(np.isfinite(chi))
    assert np.all(kappa >= 0)
    assert np.all(sigma >= 0)
    assert np.allclose(kappa + sigma, chi, rtol=2e-14, atol=0)
    assert np.allclose(sigma, scr, rtol=2e-14, atol=0)
    assert np.allclose(chi, cop, rtol=2e-14, atol=0)

    conwl5 = np.exp(50.7649141 - 5 * np.log(linelist.wlcent[0]))
    hnuk = 1.43868e8 / linelist.wlcent[0]
    planck = conwl5 / (np.exp(hnuk / atmo.temp) - 1)
    assert np.allclose(csf, planck, rtol=2e-14, atol=0)

    jbar, scattering_source = libsme.GetContinuumScatteringSource(
        linelist.wlcent[0]
    )
    assert jbar.shape == scattering_source.shape == chi.shape
    assert np.all(np.isfinite(jbar))
    assert np.all(np.isfinite(scattering_source))
    assert np.all(scattering_source > 0)
    expected_source = (kappa * planck + sigma * jbar) / chi
    assert np.allclose(scattering_source, expected_source, rtol=2e-14, atol=0)

    libsme.SetContinuumScatteringSourceMode(True)
    _, _, _, _, csf_scattering = libsme.GetLineOpacity(linelist.wlcent[0])
    assert np.allclose(csf_scattering, scattering_source, rtol=2e-14, atol=0)
    nw_scattering, wave_scattering, synth_scattering, cont_scattering = libsme.Transf(
        mu, accrt=accrt, accwi=accwt
    )
    assert nw_scattering == nw
    assert np.allclose(wave_scattering, wave)
    assert np.all(np.isfinite(synth_scattering))
    assert np.all(np.isfinite(cont_scattering))
    if np.any(sigma > 0):
        assert not np.allclose(cont_scattering, cont, rtol=1e-8, atol=0)
    libsme.SetContinuumScatteringSourceMode(False)

    libsme.GetLineRange()
    for switch in [
        "COPSTD",
        "COPRED",
        "COPBLU",
        "AHYD",
        "AH2P",
        "AHMIN",
        "SIGH",
        "AHE1",
        "AHE2",
        "AHEMIN",
        "SIGHE",
        # "ACOOL",
        # "ALUKE",
        "AHOT",
        "SIGEL",
        "SIGH2",
    ]:
        libsme.GetOpacity(switch)


@pytest.mark.parametrize("spherical", [False, True])
def test_fixed_grid_interval_index_preserves_exact_output(
    monkeypatch,
    spherical,
    libsme,
    linelist,
    teff,
    grav,
    vturb,
    atmo,
    abund,
    vw_scale,
    wfirst,
    wlast,
    mu,
):
    """Interval lookup changes line discovery, not the synthesized values."""
    model_atmo = atmo
    if spherical:
        model_atmo = deepcopy(atmo)
        model_atmo.geom = "SPH"
        model_atmo.radius = 10.0
        model_atmo.height = np.linspace(4e7, 0.0, len(model_atmo.rhox))

    libsme.SetLibraryPath()
    libsme.InputLineList(linelist)
    libsme.InputModel(teff, grav, vturb, model_atmo)
    libsme.InputAbund(abund)
    libsme.Ionization(0)
    libsme.SetVWscale(vw_scale)
    libsme.SetH2broad()
    libsme.InputWaveRange(wfirst, wlast)
    libsme.Opacity()

    threshold = 1e-4
    almax, line_range = libsme.ALMAXRange(accrt=threshold)
    strong = np.asarray(almax) >= threshold
    libsme.InputLinePrecomputedInfo(
        line_range[:, 0], line_range[:, 1], strong
    )
    libsme.SetLineInfoMode(2)
    fixed_wave = np.linspace(wfirst, wlast, 2001)

    try:
        # The first transfer consumes the one-shot LINEOP/Voigt state left by
        # ALMAXRange.  The second transfer recomputes it normally.  Equality
        # therefore checks the reuse path against the established path while
        # the interval-index setting is changed at the same time.
        monkeypatch.setenv("SME_INTERVAL_INDEX", "0")
        nw_full, wave_full, synth_full, cont_full = libsme.Transf(
            mu, wave=fixed_wave, accrt=threshold, accwi=3e-3
        )

        monkeypatch.setenv("SME_INTERVAL_INDEX", "1")
        nw_indexed, wave_indexed, synth_indexed, cont_indexed = libsme.Transf(
            mu, wave=fixed_wave, accrt=threshold, accwi=3e-3
        )

        # Any line-opacity dependency must invalidate the hand-off, even when
        # the numerical value happens to be unchanged.
        almax, line_range = libsme.ALMAXRange(accrt=threshold)
        strong = np.asarray(almax) >= threshold
        libsme.InputLinePrecomputedInfo(
            line_range[:, 0], line_range[:, 1], strong
        )
        libsme.SetVWscale(vw_scale)
        nw_invalidated, wave_invalidated, synth_invalidated, cont_invalidated = (
            libsme.Transf(mu, wave=fixed_wave, accrt=threshold, accwi=3e-3)
        )
    finally:
        libsme.SetLineInfoMode(0)

    assert nw_indexed == nw_full
    assert np.array_equal(wave_indexed, wave_full)
    assert np.array_equal(synth_indexed, synth_full)
    assert np.array_equal(cont_indexed, cont_full)
    assert nw_invalidated == nw_full
    assert np.array_equal(wave_invalidated, wave_full)
    assert np.array_equal(synth_invalidated, synth_full)
    assert np.array_equal(cont_invalidated, cont_full)


@pytest.mark.parametrize("spherical", [False, True])
def test_adaptive_grid_is_invariant_to_mu_order(
    spherical,
    libsme,
    linelist,
    teff,
    grav,
    vturb,
    atmo,
    abund,
    vw_scale,
    wfirst,
    wlast,
):
    """The adaptive-grid reference ray is the largest mu, not array index 0."""
    model_atmo = atmo
    if spherical:
        model_atmo = deepcopy(atmo)
        model_atmo.geom = "SPH"
        model_atmo.radius = 10.0
        model_atmo.height = np.linspace(4e7, 0.0, len(model_atmo.rhox))

    libsme.SetLibraryPath()
    libsme.InputLineList(linelist)
    libsme.InputModel(teff, grav, vturb, model_atmo)
    libsme.InputAbund(abund)
    libsme.Ionization(0)
    libsme.SetVWscale(vw_scale)
    libsme.SetH2broad()
    libsme.InputWaveRange(wfirst, wlast)
    libsme.Opacity()

    mu_forward = np.array([0.2, 0.6, 1.0])
    mu_reverse = mu_forward[::-1].copy()
    args = {"accrt": 1e-6, "accwi": 1e-4}

    nw_forward, wave_forward, synth_forward, cont_forward = libsme.Transf(
        mu_forward, **args
    )
    nw_reverse, wave_reverse, synth_reverse, cont_reverse = libsme.Transf(
        mu_reverse, **args
    )

    assert nw_reverse == nw_forward
    assert np.array_equal(wave_reverse, wave_forward)
    # Spherical intensities have a separate pre-existing ray-order dependence;
    # this regression scopes the spherical assertion to grid construction.
    if not spherical:
        assert np.array_equal(synth_reverse[::-1], synth_forward)
        assert np.array_equal(cont_reverse[::-1], cont_forward)


def test_continuum_scattering_source_mode_changes_spherical(
    libsme,
    linelist,
    teff,
    grav,
    vturb,
    atmo,
    abund,
    vw_scale,
    wfirst,
    wlast,
    mu,
    accrt,
    accwt,
):
    """Continuum scattering source mode works for spherical transfer."""
    spherical_atmo = deepcopy(atmo)
    spherical_atmo.geom = "SPH"
    spherical_atmo.radius = 10.0
    spherical_atmo.height = np.linspace(4e7, 0.0, len(spherical_atmo.rhox))

    libsme.SetLibraryPath()
    libsme.InputLineList(linelist)
    libsme.InputModel(teff, grav, vturb, spherical_atmo)
    libsme.InputAbund(abund)
    libsme.Ionization(0)
    libsme.SetVWscale(vw_scale)
    libsme.SetH2broad()
    libsme.InputWaveRange(wfirst, wlast)
    libsme.Opacity()

    nw, wave, synth, cont = libsme.Transf(mu, accrt=accrt, accwi=accwt)
    assert nw == len(wave) == synth.shape[-1] == cont.shape[-1]
    assert np.all(np.isfinite(synth))
    assert np.all(np.isfinite(cont))

    lop, cop, scr, tsf_legacy, csf_legacy = libsme.GetLineOpacity(
        linelist.wlcent[0]
    )
    jbar, scattering_source = libsme.GetContinuumScatteringSource(linelist.wlcent[0])
    assert jbar.shape == scattering_source.shape == cop.shape
    assert np.all(np.isfinite(jbar))
    assert np.all(np.isfinite(scattering_source))
    assert np.all(scattering_source > 0)
    assert not np.allclose(scattering_source, csf_legacy, rtol=1e-8, atol=0)

    libsme.SetContinuumScatteringSourceMode(True)
    _, _, _, tsf_scattering, csf_scattering = libsme.GetLineOpacity(
        linelist.wlcent[0]
    )
    nw_scattering, wave_scattering, synth_scattering, cont_scattering = libsme.Transf(
        mu, accrt=accrt, accwi=accwt
    )

    assert nw_scattering == nw
    assert np.allclose(wave_scattering, wave)
    assert np.all(np.isfinite(synth_scattering))
    assert np.all(np.isfinite(cont_scattering))
    assert csf_scattering.shape == scattering_source.shape
    assert np.allclose(csf_scattering, scattering_source, rtol=2e-14, atol=0)
    assert tsf_scattering.shape == tsf_legacy.shape
    assert np.all(np.isfinite(tsf_scattering))
    assert not np.allclose(cont_scattering, cont, rtol=1e-8, atol=0)
    libsme.SetContinuumScatteringSourceMode(False)
