# -*- coding: utf-8 -*-
# TODO implement NLTE tests

import os
import tempfile
from types import SimpleNamespace
from os.path import dirname

import numpy as np
import pytest

from pysme.abund import Abund
from pysme.linelist.vald import ValdFile
from pysme.nlte import DirectAccessFile, Grid, NLTE, nlte
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


def test_update_coefficients_short_format_defaults_to_warning_and_lte(caplog):
    sme = make_minimum_structure()
    sme.linelist = ValdFile("{}/testcase1.lin".format(cwd))
    sme.nlte.set_nlte("Ca", "marcs2012p_t1.0_Ca.grd")
    sme.nlte.first = True
    dll = _FakeDLL()

    with caplog.at_level("WARNING"):
        result = sme.nlte.update_coefficients(sme, dll, lfs_nlte=None)

    assert result is sme
    assert dll.reset_calls == 1
    assert "Line formation will proceed under LTE." in caplog.text


def test_update_coefficients_short_format_raises_in_strict_mode():
    sme = make_minimum_structure()
    sme.linelist = ValdFile("{}/testcase1.lin".format(cwd))
    sme.nlte.set_nlte("Ca", "marcs2012p_t1.0_Ca.grd")
    sme.nlte.strict = True
    sme.nlte.first = True
    dll = _FakeDLL()

    with pytest.raises(ValueError, match="Strict NLTE mode forbids LTE fallback"):
        sme.nlte.update_coefficients(sme, dll, lfs_nlte=None)


def test_update_coefficients_raises_in_strict_mode_when_nlte_lines_missing(monkeypatch):
    sme = make_minimum_structure()
    sme.nlte.set_nlte("Ca", "marcs2012p_t1.0_Ca.grd")
    sme.nlte.strict = True
    sme.nlte.first = True
    dll = _FakeDLL()
    fake_grid = _FakeRuntimeGrid(iused=[False, False])

    monkeypatch.setattr(NLTE, "get_grid", lambda self, sme_obj, elem, lfs: fake_grid)

    with pytest.raises(RuntimeError, match="Strict NLTE mode forbids LTE fallback"):
        sme.nlte.update_coefficients(sme, dll, lfs_nlte=None)


def test_update_coefficients_default_mode_removes_element_when_nlte_lines_missing(monkeypatch, caplog):
    sme = make_minimum_structure()
    sme.nlte.set_nlte("Ca", "marcs2012p_t1.0_Ca.grd")
    sme.nlte.first = True
    dll = _FakeDLL()
    fake_grid = _FakeRuntimeGrid(iused=[False, False])

    monkeypatch.setattr(NLTE, "get_grid", lambda self, sme_obj, elem, lfs: fake_grid)

    with caplog.at_level("WARNING"):
        sme.nlte.update_coefficients(sme, dll, lfs_nlte=None)

    assert "No Ca NLTE lines found" in caplog.text
    assert "Ca" not in sme.nlte.elements


class _FakeLineList:
    def __init__(self, species, use_indices=None):
        self.species = np.asarray(species)
        columns = []
        if use_indices is not None:
            columns.append("use_indices")
        self._lines = SimpleNamespace(columns=columns)
        self._use_indices = None if use_indices is None else np.asarray(use_indices, dtype=bool)

    def __getitem__(self, item):
        if isinstance(item, str):
            if item == "use_indices" and self._use_indices is not None:
                return self._use_indices
            raise KeyError(item)
        return _FakeLineList(self.species[item])


class _FakeSME:
    def __init__(self, linelist):
        self.linelist = linelist


class _FakeDLL:
    def __init__(self):
        self.reset_calls = 0
        self.departure_calls = []

    def ResetDepartureCoefficients(self):
        self.reset_calls += 1

    def InputDepartureCoefficients(self, bmat, lineindex):
        self.departure_calls.append((np.asarray(bmat), int(lineindex)))


class _FakeRuntimeGrid:
    def __init__(self, iused, linerefs=None, lineindices=None, bmat=None):
        self.iused = np.asarray(iused, dtype=bool)
        self.linerefs = (
            np.asarray(linerefs, dtype=int)
            if linerefs is not None
            else np.zeros((0, 2), dtype=int)
        )
        self.lineindices = (
            np.asarray(lineindices, dtype=int)
            if lineindices is not None
            else np.zeros(len(self.linerefs), dtype=int)
        )
        self._bmat = bmat

    def get(self, abund, teff, logg, monh, atmo):
        return self._bmat


def test_update_coefficients_remaps_python_indices_after_line_filtering(monkeypatch):
    sme = make_minimum_structure()
    sme.nlte.set_nlte("Ca", "marcs2012p_t1.0_Ca.grd")
    dll = _FakeDLL()
    fake_grid = _FakeRuntimeGrid(
        iused=[True, True],
        linerefs=[[0, 1], [0, 1], [0, 1]],
        lineindices=[0, 2, 3],
        bmat=np.ones((4, 2)),
    )
    monkeypatch.setattr(NLTE, "get_grid", lambda self, sme_obj, elem, lfs: fake_grid)

    # Python row 1 was discarded, so rows 0, 2 and 3 become SMElib rows 0, 1 and 2.
    sme.nlte.update_coefficients(
        sme,
        dll,
        lfs_nlte=None,
        line_index_map=np.array([0, -1, 1, 2]),
    )

    assert [lineindex for _, lineindex in dll.departure_calls] == [0, 1, 2]


def test_update_coefficients_skips_nlte_transition_discarded_by_smelib(monkeypatch):
    sme = make_minimum_structure()
    sme.nlte.set_nlte("Ca", "marcs2012p_t1.0_Ca.grd")
    dll = _FakeDLL()
    fake_grid = _FakeRuntimeGrid(
        iused=[True, True],
        linerefs=[[0, 1]],
        lineindices=[1],
        bmat=np.ones((4, 2)),
    )
    monkeypatch.setattr(NLTE, "get_grid", lambda self, sme_obj, elem, lfs: fake_grid)

    sme.nlte.update_coefficients(
        sme,
        dll,
        lfs_nlte=None,
        line_index_map=np.array([0, -1, 1]),
    )

    assert dll.departure_calls == []


def test_update_coefficients_rejects_invalid_line_index_mapping(monkeypatch):
    sme = make_minimum_structure()
    sme.nlte.set_nlte("Ca", "marcs2012p_t1.0_Ca.grd")
    dll = _FakeDLL()
    fake_grid = _FakeRuntimeGrid(
        iused=[True, True],
        linerefs=[[0, 1]],
        lineindices=[3],
        bmat=np.ones((4, 2)),
    )
    monkeypatch.setattr(NLTE, "get_grid", lambda self, sme_obj, elem, lfs: fake_grid)

    with pytest.raises(RuntimeError, match="outside the Python-to-SMElib"):
        sme.nlte.update_coefficients(
            sme,
            dll,
            lfs_nlte=None,
            line_index_map=np.array([0, -1, 1]),
        )

    with pytest.raises(ValueError, match="contiguous SMElib indices"):
        sme.nlte.update_coefficients(
            sme,
            dll,
            lfs_nlte=None,
            line_index_map=np.array([0, -1, 2, 3]),
        )

    with pytest.raises(ValueError, match="only use -1"):
        sme.nlte.update_coefficients(
            sme,
            dll,
            lfs_nlte=None,
            line_index_map=np.array([0, -2, 1]),
        )


def _make_grid_for_matching_cache(selection="energy"):
    grid = Grid.__new__(Grid)
    grid.elem = "Ti"
    grid.selection = selection
    grid.grid_name = "fake.grd"
    grid._grid_conf = np.array(["c1", "c2"])
    grid._grid_term = np.array(["t1", "t2"])
    grid._grid_species = np.array(["Ti 1", "Ti 1"])
    grid._grid_rotnum = np.array([1.0, 2.0])
    grid._grid_energies = np.array([0.1, 0.2])
    grid._grid_citation_info = "citation"
    grid.citation_info = "citation"
    grid._match_cache = {}
    grid.limits = {}
    grid.bgrid = None
    grid.depth = None
    grid.linerefs = None
    grid.lineindices = None
    grid.iused = None
    grid._active_match_key = None
    grid.first_warning = True
    return grid


def test_renew_linelist_reuses_matching_for_unchanged_use_indices(monkeypatch):
    species = np.array(["Ti 1", "Fe 1", "Ti 1"])
    use_indices = np.array([True, False, True])
    sme = _FakeSME(_FakeLineList(species, use_indices=use_indices))
    grid = _make_grid_for_matching_cache()

    calls = {"count": 0}

    def fake_select(self, conf, term, species, rotnum, energies):
        calls["count"] += 1
        return (
            np.array([0, 2]),
            np.array([[0, 1], [1, 0]]),
            np.array([True, True]),
        )

    monkeypatch.setattr(Grid, "select_energies", fake_select)

    grid.renew_linelist(sme)
    first_lineindices = grid.lineindices.copy()
    first_linerefs = grid.linerefs.copy()
    first_iused = grid.iused.copy()

    grid.renew_linelist(sme)

    assert calls["count"] == 1
    assert np.array_equal(grid.lineindices, first_lineindices)
    assert np.array_equal(grid.linerefs, first_linerefs)
    assert np.array_equal(grid.iused, first_iused)


def test_initialized_grid_skips_first_renew_for_same_use_indices(monkeypatch):
    species = np.array(["Ti 1", "Fe 1", "Ti 1"])
    use_indices = np.array([True, False, True])
    sme = _FakeSME(_FakeLineList(species, use_indices=use_indices))
    grid = _make_grid_for_matching_cache()

    calls = {"count": 0}

    def fake_select(self, conf, term, species, rotnum, energies):
        calls["count"] += 1
        return (
            np.array([0, 2]),
            np.array([[0, 1], [1, 0]]),
            np.array([True, True]),
        )

    monkeypatch.setattr(Grid, "select_energies", fake_select)

    key = grid._matching_cache_key(sme)
    grid._refresh_matching(sme, key=key)
    grid.renew_linelist(sme)

    assert calls["count"] == 1


def test_renew_linelist_recomputes_matching_when_use_indices_change(monkeypatch):
    species = np.array(["Ti 1", "Fe 1", "Ti 1"])
    grid = _make_grid_for_matching_cache()

    calls = {"count": 0}

    def fake_select(self, conf, term, species, rotnum, energies):
        calls["count"] += 1
        return (
            np.array([0]),
            np.array([[0, 0]]),
            np.array([True, False]),
        )

    monkeypatch.setattr(Grid, "select_energies", fake_select)

    sme1 = _FakeSME(_FakeLineList(species, use_indices=np.array([True, False, False])))
    sme2 = _FakeSME(_FakeLineList(species, use_indices=np.array([False, False, True])))

    grid.renew_linelist(sme1)
    grid.renew_linelist(sme2)

    assert calls["count"] == 2


def test_metal_scaled_rel_abund_is_unchanged():
    grid = _make_grid_for_abundance_test("Mg")
    abund = Abund(pattern="asplund2021", monh=0)

    expected = grid.solar_rel_abund(abund, "Mg") - grid.solar_rel_abund(abund, "Fe")

    assert grid.scaled_rel_abund(abund) == pytest.approx(expected)


# Regression tests for NLTE grid interpolation at and near grid edges (PR #39).
# A synthetic grid with a non-uniform Teff axis; departure coefficients depend only
# on Teff and are piecewise linear between the nodes, so the expected values are exact.
_EDGE_TEFF = np.array([3800.0, 3900.0, 4000.0, 4250.0])
_EDGE_B = np.array([1.0, 1.2, 1.3, 2.0])
_EDGE_DEPTH = 10.0 ** np.linspace(-4.0, 1.0, 6)


def _make_grid_for_edge_test(subgrid_teff=3):
    nlevel = 2
    teff, b = _EDGE_TEFF, _EDGE_B
    grid = Grid.__new__(Grid)
    grid.elem = "Fe"
    grid.first_warning = True
    grid.depth_name = "rhox"
    grid.iused = np.arange(nlevel)
    grid.subgrid_size = [1, subgrid_teff, 1, 1]
    grid._xfe = np.array([0.0])
    grid._teff = teff
    grid._grav = np.array([4.0])
    grid._feh = np.array([0.0])
    # keys are indexed [feh, grav, teff, abund]; each model is (level, depth)
    grid._keys = np.array([[[[f"m{i}"] for i in range(teff.size)]]])
    grid.directory = {
        f"m{i}": np.tile([[1.0], [2.0]], (1, _EDGE_DEPTH.size)) * b[i]
        for i in range(teff.size)
    }
    grid._depth = np.broadcast_to(_EDGE_DEPTH, (1, 1, teff.size, _EDGE_DEPTH.size)).copy()
    return grid


def _edge_coefficients(teff, subgrid_teff=3):
    grid = _make_grid_for_edge_test(subgrid_teff)
    grid.read_grid(0.0, teff, 4.0, 0.0)
    subgrid = grid.interpolate(0.0, teff, 4.0, 0.0, {"rhox": _EDGE_DEPTH})
    return grid, np.asarray(subgrid)


@pytest.mark.parametrize(
    "teff, expected",
    [
        (3850.0, 1.1),  # uniform interior interval, unchanged by the fix
        (4010.0, 1.3 + 0.7 * 10 / 250),  # non-uniform interval: bracket is 4000/4250
        (3700.0, 0.8),  # below the lowest node: extrapolate from 3800/3900
        (4400.0, 2.0 + 0.7 * 150 / 250),  # above the highest node: extrapolate from 4000/4250
    ],
)
def test_nlte_grid_interpolation_near_edges(teff, expected):
    grid, subgrid = _edge_coefficients(teff)

    assert subgrid.shape == (_EDGE_DEPTH.size, 2)
    np.testing.assert_allclose(subgrid[:, 0], expected, rtol=1e-12)
    np.testing.assert_allclose(subgrid[:, 1], 2 * expected, rtol=1e-12)


def test_nlte_grid_cache_includes_bracketing_nodes():
    # The three nearest nodes to 4010 K are 3800/3900/4000; 4250 must be added
    grid = _make_grid_for_edge_test(subgrid_teff=3)
    grid.read_grid(0.0, 4010.0, 4.0, 0.0)

    assert 4000.0 in grid._points[1]
    assert 4250.0 in grid._points[1]


def test_nlte_grid_warns_only_outside_grid(caplog):
    with caplog.at_level("WARNING"):
        _edge_coefficients(4010.0)
    assert "Extrapolate on Fe NLTE grid" not in caplog.text

    with caplog.at_level("WARNING"):
        _edge_coefficients(4400.0)
    assert "Extrapolate on Fe NLTE grid" in caplog.text
