# -*- coding: utf-8 -*-
"""Hydrogen-line options passed to SMElib (>= 6.13.22): h_line_profile,
h_occupation, h_occupation_neutral_scale, and the wavelength medium taken from
the line list."""
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from pysme import sme_synth
from pysme.abund import Abund
from pysme.linelist.vald import ValdFile
from pysme.sme import SME_Structure
from pysme.synthesize import (
    _compute_linelist_hash,
    _hydrogen_line_state,
    _wavelength_medium,
    synthesize_spectrum,
)

LINELIST = Path(__file__).with_name("brackett_br10_br11.lin")
BR11 = 16806.538


def test_defaults():
    sme = SME_Structure()
    assert sme.h_line_profile == "convolved"
    assert sme.h_occupation is False
    assert sme.h_occupation_neutral_scale == 0.0
    with pytest.raises(ValueError):
        sme.h_line_profile = "foo"


def test_wavelength_medium_from_the_line_list():
    assert _wavelength_medium(SimpleNamespace(linelist=SimpleNamespace(medium="air"))) == "air"
    assert _wavelength_medium(SimpleNamespace(linelist=SimpleNamespace(medium="vac"))) == "vac"
    # unknown medium: air, as for VALD
    assert _wavelength_medium(SimpleNamespace(linelist=SimpleNamespace(medium=None))) == "air"


def test_cache_identity_follows_the_hydrogen_settings():
    sme = SimpleNamespace(h_line_profile="convolved", h_occupation=False,
                          h_occupation_neutral_scale=0.0,
                          linelist=SimpleNamespace(medium="air"))
    default = _hydrogen_line_state(sme)
    sme.h_occupation = True
    occupation = _hydrogen_line_state(sme)
    sme.h_line_profile = "legacy"
    legacy = _hydrogen_line_state(sme)
    assert len({default, occupation, legacy}) == 3
    ll = ValdFile(str(LINELIST))
    assert _compute_linelist_hash(ll, default) != _compute_linelist_hash(ll, occupation)


def _br11(**options):
    sme = SME_Structure()
    sme.teff, sme.logg, sme.monh, sme.vmic, sme.vmac, sme.vsini = 5772.0, 4.44, 0.0, 1.0, 2.5, 1.8
    sme.iptype, sme.ipres = "gauss", 300000.0
    sme.atmo.method, sme.atmo.source = "grid", "marcs2012.sav"
    sme.abund = Abund(pattern="asplund2009", monh=0.0)
    sme.linelist = ValdFile(str(LINELIST))
    wave = np.arange(BR11 - 5.0, BR11 + 5.01, 0.02)
    sme.wave = [wave]
    sme.vrad_flag = sme.cscale_flag = "none"
    for key, value in options.items():
        setattr(sme, key, value)
    try:
        out = synthesize_spectrum(sme, linelist_mode="all")
    except FileNotFoundError as exc:
        pytest.skip(f"atmosphere grid unavailable: {exc}")
    flux = np.asarray(out.synth[0], dtype=float)
    return float(np.trapezoid(1.0 - flux, wave)) if hasattr(np, "trapezoid") else float(np.trapz(1.0 - flux, wave))


def test_slices_keep_the_hydrogen_line_state():
    """A sliced line list keeps the hydrogen-line state of its stored line
    information (without it the information counts as stale)."""
    ll = ValdFile(str(LINELIST))
    ll.cdr_paras = np.array([5772.0, 4.44, 0.0, 1.0])
    ll.cdr_paras_h_line_state = "profile=convolved;medium=air"
    for part in (ll[0:1], ll[[]]):
        assert part.cdr_paras_h_line_state == "profile=convolved;medium=air"


def test_options_reach_smelib_and_change_the_lines():
    """The settings are set in SMElib for the synthesis; occupation
    probabilities weaken Br11 (w(11) < 1 in the deeper layers)."""
    ew_default = _br11()
    dll = sme_synth.SME_DLL()
    assert dll.GetHydrogenOccupationMode() == 0
    assert dll.GetHydrogenProfileMode() == 1
    assert dll.GetWavelengthMedium() == "air"          # the VALD line list is in air
    ew_occ = _br11(h_occupation=True)
    assert dll.GetHydrogenOccupationMode() == 1
    assert ew_occ < ew_default
    _br11(h_line_profile="legacy")
    assert dll.GetHydrogenProfileMode() == 0
    # back to the defaults for the following tests
    _br11()
    assert dll.GetHydrogenOccupationMode() == 0 and dll.GetHydrogenProfileMode() == 1


def test_older_libraries(monkeypatch):
    """Without the SMElib functions (SMElib < 6.13.22) the defaults are accepted
    silently and other settings raise."""
    class Old:
        pass
    monkeypatch.setattr(sme_synth, "_smelib", Old())
    dll = object.__new__(sme_synth.SME_DLL)
    dll.SetHydrogenProfileMode("convolved")
    dll.SetHydrogenOccupation(False)
    dll.SetWavelengthMedium("air")
    with pytest.raises(RuntimeError):
        dll.SetHydrogenProfileMode("legacy")
    with pytest.raises(RuntimeError):
        dll.SetHydrogenOccupation(True)
