#!/usr/bin/env python
from __future__ import annotations

import argparse
import json
import subprocess
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


TEMPLATE_DIR = Path(__file__).resolve().parent / "data" / "templates"
REGRESSION_DATA = Path(__file__).resolve().parent / "data"
HALPHA_LINELIST = REGRESSION_DATA / "halpha_window_cdr_union.lin"
STELLAR_LINELIST = ROOT / "test" / "extract_stellar.lin"
DELTA_LAMBDA = 0.02


def repo_relpath(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


WINDOWS = {
    "sun_halpha": {
        "wave_range": (6552.8, 6572.8),
        "teff": 5771.0,
        "logg": 4.44,
        "monh": 0.0,
        "vmic": 1.0,
        "vmac": 0.0,
        "vsini": 0.0,
        "iptype": "gauss",
        "ipres": 47000.0,
        "linelist_path": repo_relpath(HALPHA_LINELIST),
        "nlte_elements": ["H"],
        "template_name": "sun_halpha_ref.npz",
    },
    "sun_ca5002": {
        "wave_range": (4999.8, 5003.8),
        "teff": 5771.0,
        "logg": 4.44,
        "monh": 0.0,
        "vmic": 1.0,
        "vmac": 4.19,
        "vsini": 1.6,
        "linelist_path": repo_relpath(STELLAR_LINELIST),
        "nlte_elements": ["Ca"],
        "template_name": "sun_ca5002_ref.npz",
    },
    "arcturus_halpha": {
        "wave_range": (6552.8, 6572.8),
        "teff": 4277.0,
        "logg": 1.58,
        "monh": -0.55,
        "vmic": 1.43,
        "vmac": 5.12,
        "vsini": 1.6,
        "linelist_path": repo_relpath(HALPHA_LINELIST),
        "nlte_elements": ["H"],
        "template_name": "arcturus_halpha_ref.npz",
    },
}


def get_git_rev(path: Path) -> str:
    try:
        return (
            subprocess.check_output(["git", "-C", str(path), "rev-parse", "HEAD"], text=True)
            .strip()
        )
    except Exception:
        return "unknown"


def load_linelist(cfg: dict):
    w0, w1 = cfg["wave_range"]
    ll = ValdFile(ROOT / cfg["linelist_path"])
    wl = np.asarray(ll["wlcent"], dtype=float)
    return ll[(wl >= w0 - 3.0) & (wl <= w1 + 3.0)]


def synthesize_window(cfg: dict) -> tuple[np.ndarray, np.ndarray]:
    w0, w1 = cfg["wave_range"]
    wave = np.arange(w0, w1, DELTA_LAMBDA)
    sme = SME_Structure()
    sme.teff = cfg["teff"]
    sme.logg = cfg["logg"]
    sme.monh = cfg["monh"]
    sme.vmic = cfg["vmic"]
    sme.vmac = cfg["vmac"]
    sme.vsini = cfg["vsini"]
    sme.abund = Abund(pattern="solar", monh=cfg["monh"])
    if "iptype" in cfg:
        sme.iptype = cfg["iptype"]
    if "ipres" in cfg:
        sme.ipres = cfg["ipres"]
    sme.atmo.method = "grid"
    sme.linelist = load_linelist(cfg)
    sme.wave = [wave]
    sme.normalize_by_continuum = True
    for elem in cfg["nlte_elements"]:
        sme.nlte.set_nlte(elem)
    result = synthesize_spectrum(sme)
    return np.asarray(result.wave[0], dtype=float), np.asarray(result.synth[0], dtype=float)


def write_template(name: str, cfg: dict) -> Path:
    wave, flux = synthesize_window(cfg)
    out = TEMPLATE_DIR / cfg["template_name"]
    metadata = {
        "window_key": name,
        "wave_range": list(cfg["wave_range"]),
        "deltalambda": DELTA_LAMBDA,
        "teff": cfg["teff"],
        "logg": cfg["logg"],
        "monh": cfg["monh"],
        "vmic": cfg["vmic"],
        "vmac": cfg["vmac"],
        "vsini": cfg["vsini"],
        "linelist_path": cfg["linelist_path"],
        "nlte_elements": cfg["nlte_elements"],
        "atmo_source": "default",
        "pysme_commit": get_git_rev(ROOT),
        "smelib_commit": get_git_rev(ROOT / "smelib"),
    }
    np.savez(
        out,
        wave=wave,
        flux=flux,
        metadata=json.dumps(metadata, sort_keys=True),
    )
    return out


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--window",
        choices=sorted(WINDOWS),
        default=None,
        help="Only generate a single template window.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    TEMPLATE_DIR.mkdir(parents=True, exist_ok=True)
    names = [args.window] if args.window else list(WINDOWS)
    for name in names:
        out = write_template(name, WINDOWS[name])
        print(out)


if __name__ == "__main__":
    main()
