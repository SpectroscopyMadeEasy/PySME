# -*- coding: utf-8 -*-
# TODO implement NLTE tests

import os
import tempfile
from os.path import dirname

import numpy as np
import pytest

from pysme.abund import Abund
from pysme.linelist.vald import ValdFile
from pysme.nlte import DirectAccessFile, Grid, nlte
from pysme.sme import SME_Structure as SME_Struct
from pysme.sme_synth import SME_DLL
from pysme.synthesize import Synthesizer, synthesize_spectrum

from .test_largefilestorage import lfs_atmo, lfs_nlte, skipif_lfs

cwd = dirname(__file__)


def make_minimum_structure():
    sme = SME_Struct()
    sme.teff = 5000
    sme.logg = 4.4
    sme.vmic = 1
    sme.vmac = 1
    sme.vsini = 1
    sme.abund = Abund.solar()
    sme.linelist = ValdFile("{}/testcase3.lin".format(cwd))
    sme.atmo.source = "marcs2012p_t2.0.sav"
    sme.atmo.method = "grid"

    sme.wran = [[6436, 6444]]

    return sme


def test_activate_nlte():
    sme = make_minimum_structure()

    # Make sure nothing is set yet
    assert len(sme.nlte.elements) == 0

    # Add an element
    sme.nlte.set_nlte("Ca", "marcs2012p_t1.0_Ca.grd")
    assert len(sme.nlte.elements) == 1
    assert "Ca" in sme.nlte.elements

    # Add it again, shouldn't change anything
    sme.nlte.set_nlte("Ca", "marcs2012p_t1.0_Ca.grd")
    assert len(sme.nlte.elements) == 1
    assert "Ca" in sme.nlte.elements

    # Try to remove something else
    sme.nlte.remove_nlte("Na")
    assert len(sme.nlte.elements) == 1
    assert "Ca" in sme.nlte.elements

    # Remove the original element
    sme.nlte.remove_nlte("Ca")
    assert len(sme.nlte.elements) == 0

    # Add a element with a custom grid
    sme.nlte.set_nlte("Na", "test_grid.grd")
    assert len(sme.nlte.elements) == 1
    assert "Na" in sme.nlte.elements
    assert sme.nlte.grids["Na"] == "test_grid.grd"

    # Update custom grid
    sme.nlte.set_nlte("Na", "test_grid2.grd")
    assert len(sme.nlte.elements) == 1
    assert sme.nlte.grids["Na"] == "test_grid2.grd"

    # Add element without default grid
    with pytest.raises(ValueError):
        sme.nlte.set_nlte("U")

    # with a grid it should work
    sme.nlte.set_nlte("U", "test_grid.grd")
    assert sme.nlte.grids["U"] == "test_grid.grd"


@skipif_lfs
def test_run_with_nlte():
    # NOTE sme structure must have long format for NLTE
    sme = make_minimum_structure()
    sme.nlte.set_nlte("Ca", "marcs2012p_t1.0_Ca.grd")

    sme2 = synthesize_spectrum(sme)

    assert isinstance(sme2.nlte.flags, np.ndarray)
    assert np.issubdtype(sme2.nlte.flags.dtype, np.dtype("bool"))
    assert len(sme2.nlte.flags) == len(sme2.linelist)
    assert np.any(sme2.nlte.flags)


@skipif_lfs
@pytest.mark.usefixtures("lfs_nlte")
def test_short_format_vald_raises_for_nlte(lfs_nlte):
    sme = make_minimum_structure()
    sme.linelist = ValdFile("{}/testcase1.lin".format(cwd))
    sme.nlte.set_nlte("Ca", "marcs2012p_t1.0_Ca.grd")

    with pytest.raises(ValueError, match="Short-format VALD linelists are not supported for NLTE"):
        sme.nlte.get_grid(sme, "Ca", lfs_nlte)


@skipif_lfs
@pytest.mark.usefixtures("lfs_atmo", "lfs_nlte")
def test_dll(lfs_atmo, lfs_nlte):
    sme = make_minimum_structure()
    elem = "Ca"
    sme.nlte.set_nlte(elem, "marcs2012p_t1.0_Ca.grd")

    libsme = SME_DLL()
    libsme.ResetDepartureCoefficients()

    syn = Synthesizer(None, lfs_atmo=lfs_atmo, lfs_nlte=lfs_nlte)
    sme = syn.get_atmosphere(sme)
    libsme.InputLineList(sme.linelist)
    libsme.InputModel(sme.teff, sme.logg, sme.vmic, sme.atmo)

    # This is essentially what update_depcoefs does, just for one element
    counter = 0
    bmat, linerefs, lineindices = nlte(sme, libsme, elem, lfs_nlte)
    for lr, li in zip(linerefs, lineindices):
        if lr[0] != -1 and lr[1] != -1:
            counter += 1
            libsme.InputDepartureCoefficients(bmat[:, lr], li)

    flags = libsme.GetNLTEflags()
    assert np.any(flags)
    assert np.count_nonzero(flags) == counter
    assert len(flags) == len(sme.linelist)

    idx = np.where(flags)[0][0]
    coeffs = libsme.GetDepartureCoefficients(idx)
    assert coeffs is not None

    # If we reset NLTE no flags should be set
    libsme.ResetDepartureCoefficients()
    flags = libsme.GetNLTEflags()
    assert not np.any(flags)
    assert len(flags) == len(sme.linelist)

    with pytest.raises(TypeError):
        libsme.InputDepartureCoefficients(None, 0)

    with pytest.raises(TypeError):
        libsme.InputDepartureCoefficients(bmat[:, [0, 1]], 0.1)

    with pytest.raises(ValueError):
        libsme.InputDepartureCoefficients([0, 1], 10)

    with pytest.raises(RuntimeError):
        libsme.InputDepartureCoefficients(bmat[:, [0, 1]], -10)


@pytest.fixture
def temp():
    file = tempfile.NamedTemporaryFile(delete=False)
    yield file.name
    try:
        os.remove(file)
    except:
        pass


def test_read_write_direct_access_file(temp: str):
    content = {
        "hello": "world",
        "I": ["have", "the", "high", "ground"],
        "teff": np.arange(100),
    }

    DirectAccessFile.write(temp, **content)
    daf = DirectAccessFile(temp)

    for key, value in content.items():
        vf = daf[key]
        if np.issubdtype(vf.dtype, np.dtype("S")):
            vf = np.char.decode(vf)

        assert np.all(vf == value)


def _make_grid_for_abundance_test(elem, solar_pattern="grevesse2007", abund_format="Fe=12"):
    grid = Grid.__new__(Grid)
    grid.elem = elem
    grid.abund_format = abund_format
    grid.solar = Abund(pattern=solar_pattern, monh=0)
    return grid


def test_h_scaled_rel_abund_is_zero_for_asplund2021():
    grid = _make_grid_for_abundance_test("H")
    abund = Abund(pattern="asplund2021", monh=0)

    assert grid.solar_rel_abund(abund, "H") == pytest.approx(0.0)
    assert grid.scaled_rel_abund(abund) == pytest.approx(0.0)


def test_h_scaled_rel_abund_is_zero_for_grevesse2007():
    grid = _make_grid_for_abundance_test("H")
    abund = Abund(pattern="grevesse2007", monh=0)

    assert grid.solar_rel_abund(abund, "H") == pytest.approx(0.0)
    assert grid.scaled_rel_abund(abund) == pytest.approx(0.0)


def test_h_old_scaled_rel_abund_would_show_pattern_offset():
    grid = _make_grid_for_abundance_test("H")
    abund = Abund(pattern="asplund2021", monh=0)

    old_scaled_rel_abund = grid.solar_rel_abund(abund, "H") - grid.solar_rel_abund(abund, "Fe")

    assert old_scaled_rel_abund == pytest.approx(-0.01, abs=1e-6)
    assert grid.scaled_rel_abund(abund) == pytest.approx(0.0)


def test_metal_scaled_rel_abund_is_unchanged():
    grid = _make_grid_for_abundance_test("Mg")
    abund = Abund(pattern="asplund2021", monh=0)

    expected = grid.solar_rel_abund(abund, "Mg") - grid.solar_rel_abund(abund, "Fe")

    assert grid.scaled_rel_abund(abund) == pytest.approx(expected)
