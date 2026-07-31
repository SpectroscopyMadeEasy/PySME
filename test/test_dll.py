# -*- coding: utf-8 -*-
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
    return SME_DLL()


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
    assert np.all(chi >= 0)
    assert np.allclose(kappa + sigma, chi, rtol=2e-14, atol=0)
    assert np.allclose(sigma, scr, rtol=2e-14, atol=0)
    assert np.allclose(chi, cop, rtol=2e-14, atol=0)
    jbar, scattering_source = libsme.GetContinuumScatteringSource(
        linelist.wlcent[0]
    )
    conwl5 = np.exp(50.7649141 - 5 * np.log(linelist.wlcent[0]))
    hnuk = 1.43868e8 / linelist.wlcent[0]
    planck = conwl5 / (np.exp(hnuk / atmo.temp) - 1)
    expected_source = (kappa * planck + sigma * jbar) / chi
    assert np.allclose(csf, planck, rtol=2e-14, atol=0)
    assert jbar.shape == scattering_source.shape == chi.shape
    assert np.all(np.isfinite(jbar))
    assert np.all(np.isfinite(scattering_source))
    assert np.all(scattering_source > 0)
    assert np.allclose(scattering_source, expected_source, rtol=2e-14, atol=0)
    if np.any(sigma > 0):
        assert not np.allclose(scattering_source, planck, rtol=1e-8, atol=0)
    libsme.SetContinuumScatteringSourceMode(True)
    _, _, _, tsf_scattering, csf_scattering = libsme.GetLineOpacity(
        linelist.wlcent[0]
    )
    assert np.allclose(csf_scattering, scattering_source, rtol=2e-14, atol=0)
    assert np.allclose(tsf_scattering, planck, rtol=2e-14, atol=0)
    nw_scattering, wave_scattering, synth_scattering, cont_scattering = libsme.Transf(
        mu, accrt=accrt, accwi=accwt
    )
    assert nw_scattering == nw
    assert np.allclose(wave_scattering, wave)
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
