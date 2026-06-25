from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
PYSME_SRC = ROOT / "src"
if str(PYSME_SRC) not in sys.path:
    sys.path.insert(0, str(PYSME_SRC))

from pysme.abund import Abund
from pysme.linelist.vald import ValdFile
from pysme.sme import SME_Structure
from pysme.synthesize import synthesize_spectrum
from test.conftest import skipif_smelib


pytestmark = skipif_smelib


TEMPLATE_DIR = Path(__file__).resolve().parent / "data" / "templates"
MEAN_ABS_LIMIT = 5e-4
MAX_ABS_LIMIT = 3e-3
CORE_DEPTH_LIMIT = 1e-3


def load_template(name: str) -> tuple[np.ndarray, np.ndarray, dict]:
    data = np.load(TEMPLATE_DIR / name, allow_pickle=False)
    metadata = json.loads(str(data["metadata"]))
    return data["wave"], data["flux"], metadata


def load_linelist(metadata: dict):
    w0, w1 = metadata["wave_range"]
    ll = ValdFile(ROOT / metadata["linelist_path"])
    wl = np.asarray(ll["wlcent"], dtype=float)
    return ll[(wl >= w0 - 3.0) & (wl <= w1 + 3.0)]


def synthesize_from_metadata(metadata: dict) -> np.ndarray:
    wave = np.arange(metadata["wave_range"][0], metadata["wave_range"][1], metadata["deltalambda"])
    sme = SME_Structure()
    sme.teff = metadata["teff"]
    sme.logg = metadata["logg"]
    sme.monh = metadata["monh"]
    sme.vmic = metadata["vmic"]
    sme.vmac = metadata["vmac"]
    sme.vsini = metadata["vsini"]
    sme.abund = Abund(pattern="solar", monh=metadata["monh"])
    if metadata["window_key"] == "sun_halpha":
        sme.iptype = "gauss"
        sme.ipres = 47000.0
    sme.atmo.method = "grid"
    sme.linelist = load_linelist(metadata)
    sme.wave = [wave]
    sme.normalize_by_continuum = True
    for elem in metadata["nlte_elements"]:
        sme.nlte.set_nlte(elem)
    out = synthesize_spectrum(sme)
    return np.asarray(out.synth[0], dtype=float)


def core_depth(flux: np.ndarray) -> float:
    return 1.0 - float(np.nanmin(flux))


def assert_template(name: str) -> None:
    wave_ref, flux_ref, metadata = load_template(name)
    flux_cur = synthesize_from_metadata(metadata)
    mean_abs = float(np.nanmean(np.abs(flux_cur - flux_ref)))
    max_abs = float(np.nanmax(np.abs(flux_cur - flux_ref)))
    depth_diff = abs(core_depth(flux_cur) - core_depth(flux_ref))
    assert mean_abs < MEAN_ABS_LIMIT
    assert max_abs < MAX_ABS_LIMIT
    assert depth_diff < CORE_DEPTH_LIMIT


def test_regression_sun_halpha():
    assert_template("sun_halpha_ref.npz")


def test_regression_sun_ca5002():
    assert_template("sun_ca5002_ref.npz")


def test_regression_arcturus_halpha():
    assert_template("arcturus_halpha_ref.npz")
