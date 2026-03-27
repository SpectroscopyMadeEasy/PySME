#!/usr/bin/env python
from __future__ import annotations

import argparse
import json
import pickle
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


OPTICAL_LINELIST = Path("/mnt/hard_disk/data/vald_linelists/3800_9000_long_HFS_custom_all.vlist")
H_TEMP_LINE_PKL = Path("/home/mingjie/researches/4GP/workspace/test/H_temp_line.pkl")
ARCTURUS_OPTICAL = Path(
    "/home/mingjie/software/my-py-packages/test/spectrum_stability_input/327245_melchiors_spectrum.fits"
)
TEMPLATE_DIR = Path(__file__).resolve().parent / "data" / "templates"
DELTA_LAMBDA = 0.02


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
        "linelist_mode": "pickle_segment",
        "linelist_path": str(H_TEMP_LINE_PKL),
        "nlte_elements": ["H"],
        "template_name": "sun_halpha_ref.npz",
    },
    "sun_ca6162": {
        "wave_range": (6152.173, 6172.173),
        "teff": 5771.0,
        "logg": 4.44,
        "monh": 0.0,
        "vmic": 1.0,
        "vmac": 4.19,
        "vsini": 1.6,
        "linelist_mode": "vald_segment",
        "linelist_path": str(OPTICAL_LINELIST),
        "nlte_elements": ["Ca"],
        "template_name": "sun_ca6162_ref.npz",
    },
    "arcturus_halpha": {
        "wave_range": (6552.8, 6572.8),
        "teff": 4277.0,
        "logg": 1.58,
        "monh": -0.55,
        "vmic": 1.43,
        "vmac": 5.12,
        "vsini": 1.6,
        "linelist_mode": "vald_segment",
        "linelist_path": str(OPTICAL_LINELIST),
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
    mode = cfg["linelist_mode"]
    w0, w1 = cfg["wave_range"]
    if mode == "pickle_segment":
        with H_TEMP_LINE_PKL.open("rb") as fh:
            ll = pickle.load(fh)
        if not hasattr(ll, "cdr_paras"):
            ll.cdr_paras = None
        wl = np.asarray(ll["wlcent"], dtype=float)
        return ll[(wl >= w0 - 3.0) & (wl <= w1 + 3.0)]
    if mode == "vald_segment":
        ll = ValdFile(cfg["linelist_path"])
        wl = np.asarray(ll["wlcent"], dtype=float)
        return ll[(wl >= w0 - 3.0) & (wl <= w1 + 3.0)]
    raise ValueError(f"Unknown linelist mode: {mode}")


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
        "linelist_mode": cfg["linelist_mode"],
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
