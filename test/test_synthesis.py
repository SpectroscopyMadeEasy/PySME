# -*- coding: utf-8 -*-
# TODO implement synthesis tests
import os

import numpy as np
import pandas as pd
import pytest

from pysme import util
from pysme.iliffe_vector import Iliffe_vector
from pysme.sme import SME_Structure as SME_Struct
from pysme.synthesize import Synthesizer, synthesize_spectrum
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
    def __init__(self, transf_wave=None):
        self.last_wave = "unset"
        self.transf_wave = transf_wave
        self._nlines = 0
        self.last_env = {}

    def SetLibraryPath(self):
        return None

    def InputWaveRange(self, *_):
        return None

    def Opacity(self):
        return None

    def SetLineInfoMode(self, *_):
        return None

    def InputLineList(self, linelist):
        self._nlines = len(linelist)
        return np.zeros(self._nlines, dtype=bool)

    def Transf(self, mu, accrt, accwi, keep_lineop, wave=None):
        self.last_wave = wave
        self.last_env = {
            "PYSME_H_OCCPROB_MODE": os.environ.get("PYSME_H_OCCPROB_MODE"),
            "PYSME_H_OCCPROB_FORM": os.environ.get("PYSME_H_OCCPROB_FORM"),
            "PYSME_RESAMPLE_NORM_MODE": os.environ.get("PYSME_RESAMPLE_NORM_MODE"),
            "PYSME_PROFILE_NLTE_H_CORR_MODE": os.environ.get("PYSME_PROFILE_NLTE_H_CORR_MODE"),
        }
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
        return np.zeros(0, dtype=float)

    def GetLineRange(self):
        return np.zeros((0, 2), dtype=float)

    def GetNLTEflags(self):
        return np.zeros(self._nlines, dtype=bool)


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
    assert summary["correction_construction"] == "separate"
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
    assert summary["correction_construction"] is None


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
    assert summary["correction_construction"] == "separate"


def test_profile_nlte_summary_reports_ratio_correction_mode(monkeypatch):
    dll = _DummyDLL(transf_wave=np.linspace(6550.0, 6575.0, 9))
    synth = Synthesizer(dll=dll)
    sme = _minimal_sme()
    sme.wran = [[6550.0, 6575.0]]
    sme.specific_intensities_only = False
    sme.profile_nlte.enabled = True
    sme.profile_nlte.element = "H"
    _set_minimal_species_linelist(sme, ["H 1"])

    monkeypatch.setenv("PYSME_PROFILE_NLTE_H_CORR_MODE", "ratio")
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
    assert summary["applied"] is True
    assert summary["correction_construction"] == "ratio"


def test_synthesize_segment_explicit_config_overrides_env_and_restores(monkeypatch):
    dll = _DummyDLL()
    synth = Synthesizer(dll=dll)
    sme = _minimal_sme()
    sme.specific_intensities_only = False
    sme.normalize_by_continuum = False
    sme.h_line_mode = "apply"
    sme.h_line_form = "wratio"
    sme.normalize_resample_mode = "ratio"
    sme.profile_nlte.correction_construction = "ratio"

    monkeypatch.setenv("PYSME_H_OCCPROB_MODE", "off")
    monkeypatch.setenv("PYSME_H_OCCPROB_FORM", "abs_only")
    monkeypatch.setenv("PYSME_RESAMPLE_NORM_MODE", "separate")
    monkeypatch.setenv("PYSME_PROFILE_NLTE_H_CORR_MODE", "separate")

    synth.synthesize_spectrum(
        sme,
        segments=[0],
        passAtmosphere=False,
        passNLTE=False,
    )

    assert dll.last_env["PYSME_H_OCCPROB_MODE"] == "apply"
    assert dll.last_env["PYSME_H_OCCPROB_FORM"] == "wratio"
    assert dll.last_env["PYSME_RESAMPLE_NORM_MODE"] == "ratio"
    assert dll.last_env["PYSME_PROFILE_NLTE_H_CORR_MODE"] == "ratio"
    assert os.environ.get("PYSME_H_OCCPROB_MODE") == "off"
    assert os.environ.get("PYSME_H_OCCPROB_FORM") == "abs_only"
    assert os.environ.get("PYSME_RESAMPLE_NORM_MODE") == "separate"
    assert os.environ.get("PYSME_PROFILE_NLTE_H_CORR_MODE") == "separate"


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
