# -*- coding: utf-8 -*-
"""Stehle -> TB09 blending of the Balmer profiles H10-H22 (h_tb09, h_tb09_window, h_tb09_file; SMElib
SetHydrogenTB09). The TB09 table (SYNSPEC's tremblay.dat) is not distributed, so the end-to-end test uses a
synthetic table in the same format."""
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from pysme import sme_synth
from pysme.abund import Abund
from pysme.linelist.linelist import LineList
from pysme.sme import SME_Structure
from pysme.synthesize import _compute_linelist_hash, _hydrogen_line_state, synthesize_spectrum

RH = 109677.576


def balmer_air(n):
    vac = 1e8 / (109678.764 - RH / n ** 2 - 82259.105)
    s = (1e8 / vac) ** 2
    return vac / (1.0000834213 + 2406030.0 / (1.30e10 - s) + 15997.0 / (3.89e9 - s))


def test_defaults_and_checks():
    sme = SME_Structure()
    assert sme.h_tb09 is False
    assert sme.h_tb09_window == (2.0, 1.5)
    assert sme.h_tb09_file is None
    sme.h_tb09_window = [2.5, 2.0]
    assert sme.h_tb09_window == (2.5, 2.0)
    for bad in ((1.0, 1.5), (1.0, -0.5), 3.0):
        with pytest.raises(ValueError):
            sme.h_tb09_window = bad


def test_cache_identity_follows_the_tb09_setting():
    sme = SimpleNamespace(h_line_profile="convolved", h_occupation=True, h_occupation_neutral_scale=0.0,
                          h_tb09=False, h_tb09_window=(2.0, 1.5), linelist=SimpleNamespace(medium="air"))
    off = _hydrogen_line_state(sme)
    sme.h_tb09 = True
    on = _hydrogen_line_state(sme)
    sme.h_tb09_window = (2.5, 2.0)
    other = _hydrogen_line_state(sme)
    assert len({off, on, other}) == 3 and "tb09=2,1.5" in on
    ll = LineList(pd.DataFrame({"species": ["H 1"], "wlcent": [3691.0], "excit": [10.2], "gflog": [-2.0],
                                "gamrad": [8.0], "gamqst": [2.0], "gamvw": [18.0]}), lineformat="short")
    assert _compute_linelist_hash(ll, off) != _compute_linelist_hash(ll, on)


@pytest.fixture(scope="module")
def table(tmp_path_factory):
    """Synthetic TB09 table in SYNSPEC's tremblay.dat (Lemke) format: Balmer H10-H22, one block."""
    la0, dla, mp = -3.0, 0.1, 55
    la = la0 + dla * np.arange(mp)
    s = np.exp(-(10 ** la / 0.5) ** 2) + 0.05 * (1.0 + 10 ** la) ** -2.5      # narrower than Stehle
    vals = [f"{0.0:8.4f}"] + [f"{v:8.4f}" for v in np.log10(s)]
    rec = ["".join(vals[i:i + 9]) for i in range(0, len(vals), 9)]
    nups = range(10, 23)
    out = [f"{1:6d}", f"{len(nups):12d}"]
    for nup in nups:
        out += [f"{2:12d}{nup:12d}{la0:24.15f}{10.0:24.15f}", f"{np.log10(2500.0):24.15f}{dla:24.15f}{0.5:24.15f}",
                f"{np.log10(2.0):24.15f}{mp:16d}{17:12d}{7:12d}"]
    for nup in nups:
        out.append(f"nl= 2 nu={nup:2d}")
        out += rec * (17 * 7)
    path = tmp_path_factory.mktemp("tb09") / "tremblay_synthetic.dat"
    path.write_text("\n".join(out) + "\n")
    return str(path)


def _h18(**options):
    sme = SME_Structure()
    sme.teff, sme.logg, sme.monh, sme.vmic, sme.vmac, sme.vsini = 6300.0, 4.2, -2.0, 1.5, 3.0, 0.0
    sme.iptype, sme.ipres = "gauss", 60000.0
    sme.atmo.method, sme.atmo.source = "grid", "marcs2012.sav"
    sme.abund = Abund(pattern="asplund2009", monh=-2.0)
    n = 18
    sme.linelist = LineList(pd.DataFrame({"species": ["H 1"], "wlcent": [balmer_air(n)], "excit": [10.199],
                                          "gflog": [np.log10(48.0 / n ** 3)], "gamrad": [8.0], "gamqst": [2.0],
                                          "gamvw": [float(n)]}), lineformat="short", medium="air")
    wave = np.arange(balmer_air(n) - 3.0, balmer_air(n) + 3.0, 0.02)
    sme.wave = [wave]
    sme.vrad_flag = sme.cscale_flag = "none"
    sme.h_occupation = True
    for key, value in options.items():
        setattr(sme, key, value)
    try:
        out = synthesize_spectrum(sme, linelist_mode="all")
    except FileNotFoundError as exc:
        pytest.skip(f"atmosphere grid unavailable: {exc}")
    return np.asarray(out.synth[0], dtype=float)


def test_blending_reaches_smelib(table, tmp_path):
    dll = sme_synth.SME_DLL()
    if dll.GetHydrogenTB09Mode() is None:
        pytest.skip("SMElib < 6.13.23 (no SetHydrogenTB09)")
    off = _h18()
    on = _h18(h_tb09=True, h_tb09_file=table)
    assert dll.GetHydrogenTB09Mode() == 1
    back = _h18()
    assert dll.GetHydrogenTB09Mode() == 0
    assert np.max(np.abs(on - off)) > 1e-3
    np.testing.assert_array_equal(back, off)
    with pytest.raises(RuntimeError, match="not found"):
        _h18(h_tb09=True, h_tb09_file=str(tmp_path / "missing.dat"))


def test_older_libraries(monkeypatch):
    """Without SetHydrogenTB09 (SMElib < 6.13.23) off is accepted silently and on raises."""
    class Old:
        pass
    monkeypatch.setattr(sme_synth, "_smelib", Old())
    dll = object.__new__(sme_synth.SME_DLL)
    dll.SetHydrogenTB09(0)
    assert dll.GetHydrogenTB09Mode() is None
    with pytest.raises(RuntimeError):
        dll.SetHydrogenTB09(1)
