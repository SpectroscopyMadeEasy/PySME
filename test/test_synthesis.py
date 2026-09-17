# -*- coding: utf-8 -*-
# TODO implement synthesis tests
from copy import deepcopy
import os
import numpy as np
import pandas as pd
import pytest

from pysme import util
from pysme.abund import Abund
from pysme.atmosphere.krzfile import KrzFile
from pysme.iliffe_vector import Iliffe_vector
from pysme.linelist.linelist import LineList
from pysme.sme import SME_Structure as SME_Struct
from pysme.synthesize import (
    Synthesizer,
    _build_smelib_line_index_map,
    _compute_linelist_hash,
    _compute_almax_lineinfo_for_sme,
    _load_lineinfo_cache_file,
    _temporary_brackett_convolution_env,
    synthesize_spectrum,
)
from .conftest import skipif_smelib


@skipif_smelib
def test_synthesis_simple(sme_2segments):
    sme = sme_2segments
    sme2 = synthesize_spectrum(sme)

    # Check if a result is there it has the expected data type
    assert sme2.synth is not None
    assert isinstance(sme2.synth, Iliffe_vector)
    assert np.all(sme2.synth.ravel() != 0)

    assert sme.wave is not None
    assert isinstance(sme2.wave, Iliffe_vector)
    assert np.issubdtype(sme2.wave.dtype, np.floating)
    assert np.all(sme2.wave.ravel() != 0)

    assert sme.spec is None


@skipif_smelib
def test_synthesis_segment(sme_2segments):
    sme = sme_2segments
    # Out of range
    with pytest.raises(IndexError):
        synthesize_spectrum(sme, segments=[3])

    with pytest.raises(IndexError):
        synthesize_spectrum(sme, segments=[-1])

    sme2 = synthesize_spectrum(sme, segments=[0])
    assert len(sme2.synth[0]) != 0
    assert len(sme2.synth[1]) == 0
    assert np.allclose(sme2.wran, np.array([[6550.0, 6560.0], [6560.0, 6574.0]]))

    assert len(sme2.wave[0]) != 0
    assert len(sme2.wave[1]) == 0

    assert sme2.wave.shape[0] == 2
    assert sme2.wave.shape[1][1] == 0

    orig = np.copy(sme2.synth[0])
    sme2 = synthesize_spectrum(sme2, segments=[1])

    assert sme2.wave.shape[1][0] != 0
    assert sme2.wave.shape[1][1] != 0

    assert np.all(sme2.synth[0] == orig)


class _DummyDLL:
    def __init__(self, transf_wave=None, discard_mask=None, nlte_flags=None):
        self.last_wave = "unset"
        self.transf_wave = transf_wave
        self.discard_mask = discard_mask
        self.nlte_flags = nlte_flags
        self._nlines = 0
        self.continuum_scattering_source_modes = []
        self.line_updates = []

    def SetLibraryPath(self):
        return None

    def InputWaveRange(self, *_):
        return None

    def InputModel(self, *_):
        return None

    def InputAbund(self, *_):
        return None

    def Ionization(self, *_):
        return None

    def SetVWscale(self, *_):
        return None

    def SetH2broad(self, *_):
        return None

    def Opacity(self):
        return None

    def SetLineInfoMode(self, *_):
        return None

    def InputLinePrecomputedInfo(self, *_):
        return None

    def SetContinuumScatteringSourceMode(self, mode):
        self.continuum_scattering_source_modes.append(int(mode))
        return None

    def InputLineList(self, linelist):
        if self.discard_mask is None:
            discard_mask = np.zeros(len(linelist), dtype=bool)
        else:
            discard_mask = np.asarray(self.discard_mask, dtype=bool)
            if discard_mask.shape != (len(linelist),):
                raise ValueError("test discard mask does not match line list")
        self._nlines = int(np.count_nonzero(~discard_mask))
        return discard_mask

    def UpdateLineList(self, atomic, species, index):
        self.line_updates.append(
            (
                np.asarray(atomic).copy(),
                np.asarray(species).copy(),
                np.asarray(index).copy(),
            )
        )

    def Transf(self, mu, accrt, accwi, keep_lineop, wave=None):
        self.last_wave = wave
        self.brackett_mode = os.environ.get("PYSME_H_STARK_CONVOLUTION")
        if wave is None:
            if self.transf_wave is None:
                wint = np.linspace(5000.0, 5001.0, 5)
            else:
                wint = np.asarray(self.transf_wave, dtype=float)
        else:
            wint = np.asarray(wave, dtype=float)
        sint = np.ones((len(mu), len(wint)), dtype=float)
        cint = np.ones((len(mu), len(wint)), dtype=float)
        return len(wint), wint, sint, cint

    def CentralDepth(self, mu, accrt):
        self.central_depth_mode = os.environ.get("PYSME_H_STARK_CONVOLUTION")
        return np.zeros(0, dtype=float)

    def GetLineRange(self):
        return np.zeros((0, 2), dtype=float)

    def ALMAXRange(self, accrt):
        self.almax_mode = os.environ.get("PYSME_H_STARK_CONVOLUTION")
        return np.zeros(self._nlines), np.zeros((self._nlines, 2))

    def GetNLTEflags(self):
        if self.nlte_flags is None:
            return np.zeros(self._nlines, dtype=bool)
        flags = np.asarray(self.nlte_flags, dtype=bool)
        if flags.shape != (self._nlines,):
            raise ValueError("test NLTE flags do not match retained lines")
        return flags


def _minimal_sme():
    sme = SME_Struct()
    sme.wran = [[5000.0, 5001.0]]
    sme.vsini = 0.0
    sme.vmac = 0.0
    sme.vrad_flag = "none"
    sme.specific_intensities_only = True
    sme.normalize_by_continuum = True
    return sme


def _set_minimal_species_linelist(sme, species):
    sme.linelist._lines = pd.DataFrame(
        {
            "species": list(species),
            "wlcent": np.full(len(species), 6562.8, dtype=float),
        }
    )
    sme.line_ion_mask = np.zeros(len(species), dtype=bool)


def _continuum_scattering_sme(cwd, spherical=False):
    sme = SME_Struct()
    sme.teff = 5750
    sme.logg = 4.5
    sme.vmic = 2.0
    sme.vmac = 0.0
    sme.vsini = 0.0
    sme.abund = Abund(monh=0, pattern="asplund2009")
    sme.linelist = LineList()
    sme.linelist.add("Fe 1", 5502.9931, 0.9582, -3.047, 7.19, -6.22, 239.249)
    sme.linelist.add("Cr 2", 5503.5955, 4.1682, -2.117, 8.37, -6.49, 195.248)
    sme.atmo = KrzFile(os.path.join(cwd, "testatmo1.krz"))
    sme.atmo.method = "embedded"
    if spherical:
        sme.atmo.geom = "SPH"
        sme.atmo.radius = 10.0
        sme.atmo.height = np.linspace(4e7, 0.0, len(sme.atmo.rhox))
    else:
        sme.atmo.geom = "PP"
    sme.wran = [[5500.0, 5600.0]]
    sme.wint = [np.linspace(5500.0, 5600.0, 41)]
    sme.mu = [1.0]
    sme.vrad_flag = "none"
    sme.cscale_flag = "none"
    sme.specific_intensities_only = True
    return sme


def _synthesize_continuum_with_mode(synth, sme, mode):
    one = deepcopy(sme)
    one.continuum_scattering_source = mode
    out = synth.synthesize_spectrum(one, passNLTE=False)
    return np.asarray(out.wint[0]), np.asarray(out.cint[0])


def test_build_smelib_line_index_map_compacts_retained_rows():
    mapping = _build_smelib_line_index_map([False, True, False, True, False])

    assert np.array_equal(mapping, np.array([0, -1, 1, -1, 2]))


def test_synthesis_passes_compact_line_index_map_to_nlte(monkeypatch):
    dll = _DummyDLL(discard_mask=[False, True, False, False])
    synth = Synthesizer(dll=dll)
    sme = _minimal_sme()
    _set_minimal_species_linelist(sme, ["Fe 1", "Ni 7", "Li 1", "Ti 1"])
    captured = {}

    def capture_update(self, sme_obj, dll_obj, lfs_nlte, line_index_map=None):
        captured["line_index_map"] = np.asarray(line_index_map).copy()
        return sme_obj

    monkeypatch.setattr(type(sme.nlte), "update_coefficients", capture_update)

    synth.synthesize_spectrum(
        sme,
        passAtmosphere=False,
        passNLTE=True,
        updateStructure=False,
    )

    assert np.array_equal(captured["line_index_map"], np.array([0, -1, 1, 2]))
    assert np.array_equal(sme.line_ion_mask, [False, True, False, False])


def test_reused_linelist_maps_incremental_updates_to_smelib_indices():
    dll = _DummyDLL(discard_mask=[False, True, False, False])
    synth = Synthesizer(dll=dll)
    sme = _minimal_sme()
    _set_minimal_species_linelist(sme, ["Fe 1", "Ni 7", "Li 1", "Ti 1"])

    synth.synthesize_spectrum(
        sme,
        passAtmosphere=False,
        passNLTE=False,
        updateStructure=False,
    )
    synth.synthesize_spectrum(
        sme,
        passLineList=False,
        passAtmosphere=False,
        passNLTE=False,
        updateStructure=False,
        updateLineList=[1, 2, 3],
    )

    assert len(dll.line_updates) == 1
    atomic, species, indices = dll.line_updates[0]
    assert atomic.shape[0] == 3
    assert np.array_equal(species, ["Fe 1", "Li 1", "Ti 1"])
    assert np.array_equal(indices, [1, 2])


def test_dynamic_synthesis_maps_selected_and_discarded_lines(monkeypatch):
    dll = _DummyDLL(
        discard_mask=[True, False, False],
        nlte_flags=[True, False],
    )
    synth = Synthesizer(dll=dll)
    sme = _minimal_sme()
    sme.specific_intensities_only = False
    sme.wave = np.linspace(5000.0, 5001.0, 11)
    _set_minimal_species_linelist(sme, ["Ni 7", "Li 1", "Ti 1", "Fe 1"])
    sme.linelist._lines["central_depth"] = [0.2, 0.2, 0.2, 0.0]
    sme.linelist._lines["line_range_s"] = [4999.0, 4999.0, 4999.0, np.nan]
    sme.linelist._lines["line_range_e"] = [5002.0, 5002.0, 5002.0, np.nan]
    sme.linelist._lines["strong"] = [True, True, True, False]
    sme.line_select_method = "cdr"
    sme.ipres = 1e9
    sme.linelist.cdr_paras = np.array([sme.teff, sme.logg, sme.monh, sme.vmic])
    sme.linelist.cdr_paras_h_stark_convolution = "legacy"
    sme.linelist.cdr_paras_thres["strong_depth"] = sme.strong_depth_thres
    sme.linelist.cdr_paras_thres["strong_bin_width"] = sme.strong_bin_width
    captured = {}

    def capture_update(self, sme_obj, dll_obj, lfs_nlte, line_index_map=None):
        captured["line_index_map"] = np.asarray(line_index_map).copy()
        return sme_obj

    monkeypatch.setattr(type(sme.nlte), "update_coefficients", capture_update)

    synth.synthesize_spectrum(
        sme,
        passAtmosphere=False,
        passNLTE=True,
        updateStructure=True,
        linelist_mode="dynamic",
    )

    assert np.array_equal(sme.linelist["use_indices"], [True, True, True, False])
    assert np.array_equal(captured["line_index_map"], [-1, 0, 1])
    assert np.array_equal(sme.line_ion_mask, [True, False, False, False])
    assert np.array_equal(sme.nlte.flags, [False, True, False, False])
    assert np.array_equal(sme.linelist["nlte_flag"], [-1, 1, 0, -1])

    captured.clear()
    synth.synthesize_spectrum(
        sme,
        passLineList=False,
        passAtmosphere=False,
        passNLTE=True,
        updateStructure=False,
    )
    assert np.array_equal(captured["line_index_map"], [-1, 0, 1])


def test_synthesize_segment_prefers_user_wint_over_cache():
    dll = _DummyDLL()
    synth = Synthesizer(dll=dll)
    sme = _minimal_sme()

    user_wint = np.linspace(5000.0, 5001.0, 7)
    cached_wint = np.linspace(5000.0, 5001.0, 9)
    sme.wint = user_wint
    synth.wint[0] = cached_wint.copy()

    synth.synthesize_segment(sme, 0, reuse_wavelength_grid=True)

    assert np.allclose(dll.last_wave, user_wint)
    assert np.allclose(synth.wint[0], cached_wint)


def test_synthesize_segment_uses_cache_when_user_wint_missing():
    dll = _DummyDLL()
    synth = Synthesizer(dll=dll)
    sme = _minimal_sme()

    cached_wint = np.linspace(5000.0, 5001.0, 9)
    synth.wint[0] = cached_wint

    synth.synthesize_segment(sme, 0, reuse_wavelength_grid=True)

    assert np.allclose(dll.last_wave, cached_wint)


def test_synthesize_segment_populates_cache_when_no_wint_available():
    dll = _DummyDLL()
    synth = Synthesizer(dll=dll)
    sme = _minimal_sme()

    synth.synthesize_segment(sme, 0, reuse_wavelength_grid=True)

    assert dll.last_wave is None
    assert 0 in synth.wint
    assert np.allclose(synth.wint[0], np.linspace(5000.0, 5001.0, 5))


def test_brackett_mode_is_explicit_and_environment_is_restored(monkeypatch):
    dll = _DummyDLL()
    synth = Synthesizer(dll=dll)
    sme = _minimal_sme()
    sme.h_stark_convolution = "legacy"
    monkeypatch.setenv("PYSME_H_STARK_CONVOLUTION", "convolution")

    synth.synthesize_segment(sme, 0)

    assert dll.brackett_mode == "legacy"
    assert dll.central_depth_mode == "legacy"
    assert os.environ["PYSME_H_STARK_CONVOLUTION"] == "convolution"


def test_brackett_environment_is_restored_after_exception(monkeypatch):
    sme = _minimal_sme()
    sme.h_stark_convolution = "convolution"
    monkeypatch.setenv("PYSME_H_STARK_CONVOLUTION", "legacy")

    with pytest.raises(RuntimeError):
        with _temporary_brackett_convolution_env(sme):
            assert os.environ["PYSME_H_STARK_CONVOLUTION"] == "convolution"
            raise RuntimeError("native synthesis failed")

    assert os.environ["PYSME_H_STARK_CONVOLUTION"] == "legacy"


@pytest.mark.parametrize(
    "external_value,expected",
    [("convolution", "convolution"), ("convolutionXYZ", "legacy")],
)
def test_brackett_none_uses_only_exact_external_mode(monkeypatch, external_value, expected):
    sme = _minimal_sme()
    sme.h_stark_convolution = None
    monkeypatch.setenv("PYSME_H_STARK_CONVOLUTION", external_value)

    with _temporary_brackett_convolution_env(sme):
        assert os.environ["PYSME_H_STARK_CONVOLUTION"] == expected

    assert os.environ["PYSME_H_STARK_CONVOLUTION"] == external_value


def test_almax_worker_scopes_brackett_mode(monkeypatch):
    dll = _DummyDLL()
    monkeypatch.setattr(Synthesizer, "get_dll", lambda self, dll_id=None: dll)
    monkeypatch.setattr(Synthesizer, "get_atmosphere", lambda self, sme: sme)
    sme = _minimal_sme()
    sme.h_stark_convolution = "convolution"
    sme.atmo.method = "embedded"
    sme.linelist._lines = pd.DataFrame(
        {"species": ["H 1"], "wlcent": [5000.0]}
    )
    monkeypatch.setenv("PYSME_H_STARK_CONVOLUTION", "legacy")

    _compute_almax_lineinfo_for_sme(sme)

    assert dll.almax_mode == "convolution"
    assert os.environ["PYSME_H_STARK_CONVOLUTION"] == "legacy"


def test_convolution_lineinfo_cache_is_namespaced_and_mode_checked(tmp_path):
    sme = _minimal_sme()
    sme.linelist._lines = pd.DataFrame(
        {"species": ["H 1"], "wlcent": [5000.0], "atomic": [1.0]}
    )
    legacy_hash = _compute_linelist_hash(sme.linelist, "legacy")
    convolution_hash = _compute_linelist_hash(sme.linelist, "convolution")
    assert legacy_hash == _compute_linelist_hash(sme.linelist)
    assert convolution_hash != legacy_hash

    cache = tmp_path / "cdr_cache.npz"
    np.savez_compressed(
        cache,
        line_info=np.array([[0.0, 0.1, 4999.9, 5000.1]]),
        method=np.array("cdr"),
        linelist_hash=np.array(legacy_hash),
        n_lines_total=np.array(1),
        h_stark_convolution=np.array("legacy"),
    )
    _load_lineinfo_cache_file(cache, "cdr", legacy_hash, 1, "legacy")
    with pytest.raises(ValueError, match="h_stark_convolution mismatch"):
        _load_lineinfo_cache_file(cache, "cdr", legacy_hash, 1, "convolution")


@pytest.mark.parametrize("method,metric", [("cdr", "central_depth"), ("almax", "almax_ratio")])
def test_switching_brackett_mode_invalidates_lineinfo(monkeypatch, method, metric):
    sme = _minimal_sme()
    sme.h_stark_convolution = "convolution"
    sme.line_select_method = method
    sme.line_select_policy = "strict"
    sme.line_select_recompute = "if_stale"
    sme.linelist._lines = pd.DataFrame(
        {
            "species": ["H 1"],
            "wlcent": [5000.0],
            metric: [0.02],
            "line_range_s": [4999.9],
            "line_range_e": [5000.1],
            "strong": [True],
        }
    )
    if method == "cdr":
        sme.linelist.cdr_paras = np.array([sme.teff, sme.logg, sme.monh, sme.vmic])
        sme.linelist.cdr_paras_h_stark_convolution = "legacy"
    else:
        sme.linelist.almax_paras = np.array(
            [sme.teff, sme.logg, sme.monh, sme.vmic, sme.accrt, sme.accrt, 0.0, 0.2]
        )
        sme.linelist.almax_paras_h_stark_convolution = "legacy"

    called = {"count": 0}

    def fail_if_not_invalidated(self, sme, **kwargs):
        called["count"] += 1
        raise RuntimeError("mode switch invalidated line info")

    monkeypatch.setattr(Synthesizer, f"update_{method}", fail_if_not_invalidated)
    synth = Synthesizer(dll=_DummyDLL())

    with pytest.raises(RuntimeError, match="mode switch"):
        synth.synthesize_spectrum(
            sme,
            passAtmosphere=False,
            passNLTE=False,
            updateStructure=False,
        )
    assert called["count"] == 1


def test_specific_intensities_only_updates_sme_and_trims_to_wran():
    dll = _DummyDLL(transf_wave=np.linspace(4999.5, 5001.5, 9))
    synth = Synthesizer(dll=dll)
    sme = _minimal_sme()

    out = synth.synthesize_spectrum(
        sme,
        segments=[0],
        passAtmosphere=False,
        passNLTE=False,
    )

    assert out is sme
    assert hasattr(sme, "wint")
    assert hasattr(sme, "sint")
    assert hasattr(sme, "cint")

    w = np.asarray(sme.wint[0], dtype=float)
    sint = np.asarray(sme.sint[0], dtype=float)
    cint = np.asarray(sme.cint[0], dtype=float)

    assert w.size > 0
    assert sint.size > 0
    assert cint.size > 0
    assert np.all(w >= 5000.0)
    assert np.all(w <= 5001.0)
    assert w.size < 9
    assert sint.shape == (len(sme.mu), w.size)
    assert cint.shape == (len(sme.mu), w.size)


def test_synthesize_spectrum_sets_continuum_scattering_source_mode_each_time():
    dll = _DummyDLL(transf_wave=np.linspace(5000.0, 5001.0, 5))
    synth = Synthesizer(dll=dll)
    sme = _minimal_sme()

    sme.continuum_scattering_source = False
    synth.synthesize_spectrum(
        sme,
        segments=[0],
        passAtmosphere=False,
        passNLTE=False,
    )
    sme.continuum_scattering_source = True
    synth.synthesize_spectrum(
        sme,
        segments=[0],
        passAtmosphere=False,
        passNLTE=False,
    )
    sme.continuum_scattering_source = False
    synth.synthesize_spectrum(
        sme,
        segments=[0],
        passAtmosphere=False,
        passNLTE=False,
    )

    assert dll.continuum_scattering_source_modes == [0, 1, 0]


@skipif_smelib
def test_high_level_continuum_scattering_source_off_on_off_pp():
    sme = _continuum_scattering_sme(os.path.dirname(__file__), spherical=False)
    synth = Synthesizer()

    wave_off_1, cont_off_1 = _synthesize_continuum_with_mode(synth, sme, False)
    wave_on, cont_on = _synthesize_continuum_with_mode(synth, sme, True)
    wave_off_2, cont_off_2 = _synthesize_continuum_with_mode(synth, sme, False)

    assert np.allclose(wave_on, wave_off_1)
    assert np.allclose(wave_off_2, wave_off_1)
    assert np.all(np.isfinite(cont_off_1))
    assert np.all(np.isfinite(cont_on))
    assert np.all(np.isfinite(cont_off_2))
    assert np.allclose(cont_off_2, cont_off_1, rtol=0, atol=0)
    assert not np.allclose(cont_on, cont_off_1, rtol=1e-8, atol=0)


@skipif_smelib
def test_high_level_continuum_scattering_source_spherical_on():
    sme = _continuum_scattering_sme(os.path.dirname(__file__), spherical=True)
    synth = Synthesizer()

    wave_off, cont_off = _synthesize_continuum_with_mode(synth, sme, False)
    wave_on, cont_on = _synthesize_continuum_with_mode(synth, sme, True)
    wave_off_2, cont_off_2 = _synthesize_continuum_with_mode(synth, sme, False)

    assert np.allclose(wave_on, wave_off)
    assert np.allclose(wave_off_2, wave_off)
    assert cont_on.shape == cont_off.shape
    assert np.all(np.isfinite(cont_off))
    assert np.all(np.isfinite(cont_on))
    assert np.all(np.isfinite(cont_off_2))
    assert np.allclose(cont_off_2, cont_off, rtol=0, atol=0)
    assert not np.allclose(cont_on, cont_off, rtol=1e-8, atol=0)


def test_profile_nlte_h_summary_uses_default_provider(monkeypatch):
    dll = _DummyDLL(transf_wave=np.linspace(6550.0, 6575.0, 9))
    synth = Synthesizer(dll=dll)
    sme = _minimal_sme()
    sme.wran = [[6550.0, 6575.0]]
    sme.specific_intensities_only = False
    sme.profile_nlte.enabled = True
    sme.profile_nlte.element = "H"
    _set_minimal_species_linelist(sme, ["H 1"])

    monkeypatch.setattr(
        Synthesizer,
        "get_H_3dnlte_correction_rbf",
        lambda self, sme: np.ones((len(sme.mu), len(np.asarray(util.lambda_H_3DNLTE)))),
    )

    out = synth.synthesize_spectrum(
        sme,
        segments=[0],
        passAtmosphere=False,
        passNLTE=False,
    )

    summary = out.profile_nlte.summary
    assert summary["requested"] is True
    assert summary["applied"] is True
    assert summary["element"] == "H"
    assert summary["provider"] == "pysme_h_3dnlte_rbf"
    assert summary["species"] == "H 1"
    assert summary["profile_kind"] == "intensity_ratio"
    assert summary["data_key"] == "data.hlineprof"
    assert summary["data_source"] == "lineprof.dat"
    assert summary["parameter_axes"] == ["teff", "logg", "monh", "mu"]
    assert summary["reference_label"] == "Amarsi et al. (2018, A&A, 615, A139)"
    assert "Amarsi2018Balmer3DNLTE" in summary["citation_info"]
    assert summary["fallback_reason"] is None
    assert summary["applied_windows_air"] == [[6550.0, 6575.0]]


def test_profile_nlte_summary_skips_when_species_missing():
    dll = _DummyDLL(transf_wave=np.linspace(6550.0, 6575.0, 9))
    synth = Synthesizer(dll=dll)
    sme = _minimal_sme()
    sme.wran = [[6550.0, 6575.0]]
    sme.specific_intensities_only = False
    sme.profile_nlte.enabled = True
    sme.profile_nlte.element = "H"
    _set_minimal_species_linelist(sme, ["Fe 1"])

    out = synth.synthesize_spectrum(
        sme,
        segments=[0],
        passAtmosphere=False,
        passNLTE=False,
    )

    summary = out.profile_nlte.summary
    assert summary["requested"] is True
    assert summary["applied"] is False
    assert summary["fallback"] is True
    assert summary["fallback_reason"] == "no_matching_species_in_linelist"


def test_profile_nlte_warning_when_species_missing(caplog):
    dll = _DummyDLL(transf_wave=np.linspace(6550.0, 6575.0, 9))
    synth = Synthesizer(dll=dll)
    sme = _minimal_sme()
    sme.wran = [[6550.0, 6575.0]]
    sme.specific_intensities_only = False
    sme.profile_nlte.enabled = True
    sme.profile_nlte.element = "H"
    _set_minimal_species_linelist(sme, ["Fe 1"])

    with caplog.at_level("WARNING", logger="pysme.synthesize"):
        synth.synthesize_spectrum(
            sme,
            segments=[0],
            passAtmosphere=False,
            passNLTE=False,
        )

    assert "was requested but not applied" in caplog.text
    assert "does not contain species 'H 1'" in caplog.text


def test_profile_nlte_legacy_tdnlte_h_maps_to_profile_summary(monkeypatch):
    dll = _DummyDLL(transf_wave=np.linspace(6550.0, 6575.0, 9))
    synth = Synthesizer(dll=dll)
    sme = _minimal_sme()
    sme.wran = [[6550.0, 6575.0]]
    sme.specific_intensities_only = False
    sme.tdnlte_H = True
    _set_minimal_species_linelist(sme, ["H 1"])

    monkeypatch.setattr(
        Synthesizer,
        "get_H_3dnlte_correction_rbf",
        lambda self, sme: np.ones((len(sme.mu), len(np.asarray(util.lambda_H_3DNLTE)))),
    )

    out = synth.synthesize_spectrum(
        sme,
        segments=[0],
        passAtmosphere=False,
        passNLTE=False,
    )

    summary = out.profile_nlte.summary
    assert summary["requested"] is True
    assert summary["applied"] is True
    assert summary["element"] == "H"
    assert summary["provider"] == "pysme_h_3dnlte_rbf"


def _flag_strong_lines_by_bins_reference(wl, depth, bin_width=0.2, threshold=0.001):
    wl = np.asarray(wl, dtype=float)
    depth = np.asarray(depth, dtype=float)

    invalid_depth = ~np.isfinite(depth)
    depth_sanitized = np.where(invalid_depth, 0.0, depth)
    depth_sanitized = np.where(depth_sanitized > 0, depth_sanitized, 0.0)

    wl_min, wl_max = wl.min(), wl.max()
    edges = np.arange(wl_min, wl_max + bin_width, bin_width)
    bin_idx = np.searchsorted(edges, wl, side="right") - 1

    order = np.lexsort((depth_sanitized, bin_idx))
    bin_sorted = bin_idx[order]
    depth_sorted = depth_sanitized[order]

    starts = np.r_[0, np.flatnonzero(np.diff(bin_sorted)) + 1]
    ends = np.r_[starts[1:], len(depth_sorted)]

    keep_sorted = np.empty_like(depth_sorted, dtype=bool)
    for s, e in zip(starts, ends):
        if s == e:
            continue
        local = np.cumsum(depth_sorted[s:e])
        cut = np.searchsorted(local, threshold, side="right")
        keep_sorted[s:e] = True
        if cut > 0:
            keep_sorted[s : s + cut] = False

    keep = np.empty_like(keep_sorted)
    keep[order] = keep_sorted
    keep[invalid_depth] = False
    return keep


def test_flag_strong_lines_by_bins_matches_reference():
    rng = np.random.default_rng(42)
    n = 5000
    wl = rng.uniform(6000.0, 6030.0, n)
    depth = rng.uniform(-1e-4, 2e-3, n)

    # Inject special values and repeated wavelengths to stress boundary handling.
    depth[rng.choice(n, size=300, replace=False)] = 0.0
    depth[rng.choice(n, size=100, replace=False)] = np.nan
    wl[:20] = np.linspace(6005.0, 6006.9, 20)

    for bin_width in [0.2, 0.03]:
        for threshold in [1e-3, 0.0, -1e-4]:
            got = Synthesizer.flag_strong_lines_by_bins(
                wl, depth, bin_width=bin_width, threshold=threshold
            )
            ref = _flag_strong_lines_by_bins_reference(
                wl, depth, bin_width=bin_width, threshold=threshold
            )
            assert np.array_equal(got, ref)


def test_flag_strong_lines_by_bins_valid_mask_matches_sliced_call():
    rng = np.random.default_rng(7)
    n = 8000
    wl = rng.uniform(5000.0, 5100.0, n)
    depth = rng.uniform(-5e-5, 2e-3, n)
    depth[rng.choice(n, size=200, replace=False)] = np.nan
    valid_mask = rng.random(n) > 0.3

    got = Synthesizer.flag_strong_lines_by_bins(
        wl,
        depth,
        bin_width=0.2,
        threshold=1e-3,
        valid_mask=valid_mask,
    )
    ref_sub = _flag_strong_lines_by_bins_reference(
        wl[valid_mask],
        depth[valid_mask],
        bin_width=0.2,
        threshold=1e-3,
    )
    ref = np.zeros(n, dtype=bool)
    ref[valid_mask] = ref_sub

    assert np.array_equal(got, ref)
