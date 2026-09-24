"""Bounded benchmark for a possible fixed-grid line-centric opacity path.

The first phase profiles the existing wavelength-centric implementation on
10, 200, and 800 A windows.  It deliberately includes atmosphere/EOS,
ALMAX/range discovery, and transfer in the reported synthesis wall time, but
excludes parsing the large VALD file.  Native ``SME_TIMING`` output can be
enabled externally to obtain the Transf phase breakdown.

This is an analysis helper, not a public PySME interface.
"""

from __future__ import annotations

import argparse
import copy
import json
import time
from pathlib import Path

import numpy as np

from pysme.abund import Abund
from pysme.linelist.vald import ValdFile
from pysme.sme import SME_Structure
from pysme.synthesize import Synthesizer


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_SOURCE = Path(
    "/Users/mingjie/Documents/Research/data/merged_3700-9500_hfs.lin"
)
DEFAULT_CACHE = Path("/tmp/pysme_line_centric_windows")
DEFAULT_OUTPUT = ROOT / "analysis" / "line_centric_wide_benchmark.json"
CENTER_A = 5200.0
PADDING_A = 150.0
STEP_A = 0.05
WIDTHS_A = (10.0, 200.0, 800.0)
MODELS = (
    ("solar_dwarf", 5777.0, 4.44, 0.0, 1.0),
    ("cool_metal_rich_dwarf", 4250.0, 4.5, 0.3, 1.0),
)


def _window_bounds(width: float) -> tuple[float, float]:
    return CENTER_A - width / 2, CENTER_A + width / 2


def extract_vald_window(source: Path, destination: Path, lo: float, hi: float) -> int:
    """Stream one four-line-per-transition VALD window into a small file."""
    destination.parent.mkdir(parents=True, exist_ok=True)
    selected = 0
    footer: list[str] = []
    body = destination.with_suffix(destination.suffix + ".body")
    with source.open("r", encoding="utf-8", errors="replace") as src:
        next(src)
        header_2 = next(src)
        header_3 = next(src)
        with body.open("w", encoding="utf-8") as dst:
            while True:
                record = [src.readline() for _ in range(4)]
                if not record[0]:
                    break
                if any(item == "" for item in record):
                    raise RuntimeError(f"Truncated VALD record in {source}")
                fields = record[0].split(",", 2)
                if len(fields) < 3:
                    footer = record + list(src)
                    break
                try:
                    wavelength = float(fields[1])
                except ValueError:
                    # VALD files append an abundance/reference footer after
                    # the transition blocks (normally starting "undefined").
                    footer = record + list(src)
                    break
                if lo <= wavelength <= hi:
                    dst.writelines(record)
                    selected += 1

    first = (
        f"{lo:11.4f},{hi:11.4f},{selected:d},{selected:d}, 0.0, "
        "Wavelength region, Lines selected, Lines processed, Vmicro\n"
    )
    with destination.open("w", encoding="utf-8") as dst, body.open(
        "r", encoding="utf-8"
    ) as src:
        dst.write(first)
        dst.write(header_2)
        dst.write(header_3)
        for chunk in iter(lambda: src.read(1024 * 1024), ""):
            dst.write(chunk)
        dst.writelines(footer)
    body.unlink()
    return selected


def ensure_window(source: Path, cache: Path, width: float) -> Path:
    lo, hi = _window_bounds(width)
    extract_lo, extract_hi = lo - PADDING_A, hi + PADDING_A
    path = cache / f"window_{int(width):04d}A_pad_{int(PADDING_A):03d}A.lin"
    if not path.exists():
        started = time.perf_counter()
        count = extract_vald_window(source, path, extract_lo, extract_hi)
        print(
            f"extracted width={width:g} A: {count:,} lines in "
            f"{time.perf_counter() - started:.2f} s",
            flush=True,
        )
    return path


def make_sme(model: tuple, linelist: ValdFile, width: float) -> SME_Structure:
    _, teff, logg, monh, vmic = model
    lo, hi = _window_bounds(width)
    wave = np.linspace(lo, hi, int(round(width / STEP_A)) + 1)
    sme = SME_Structure()
    sme.teff = teff
    sme.logg = logg
    sme.monh = monh
    sme.vmic = vmic
    sme.vmac = 0.0
    sme.vsini = 0.0
    sme.ipres = 1.0e7
    sme.abund = Abund(monh=monh, pattern="asplund2009")
    sme.linelist = copy.deepcopy(linelist)
    sme.atmo.source = "marcs2012.sav"
    sme.atmo.method = "grid"
    sme.wave = [wave]
    sme.wint = [wave]
    sme.vrad_flag = "none"
    sme.cscale_flag = "none"
    sme.accrt = 1e-4
    sme.line_select_method = "almax"
    sme.line_select_policy = "auto"
    sme.line_select_parallel = False
    sme.line_select_recompute = "always"
    sme.continuum_grid = "adaptive"
    sme.continuum_grid_base_step = 1.0
    sme.continuum_grid_rtol = 1e-3
    sme.continuum_grid_min_step = 1e-3
    return sme


def run(args: argparse.Namespace) -> dict:
    output = {
        "center_A": CENTER_A,
        "step_A": STEP_A,
        "padding_A": PADDING_A,
        "source": str(args.source),
        "models": {},
    }
    for width in args.widths:
        path = ensure_window(args.source, args.cache, width)
        started = time.perf_counter()
        linelist = ValdFile(str(path))
        parse_sec = time.perf_counter() - started
        print(
            f"loaded width={width:g} A: {len(linelist):,} lines in {parse_sec:.2f} s",
            flush=True,
        )
        for model in MODELS:
            name = model[0]
            sme = make_sme(model, linelist, width)
            started = time.perf_counter()
            Synthesizer().synthesize_spectrum(
                sme,
                linelist_mode="all",
                smelib_lineinfo_mode=0,
            )
            wall = time.perf_counter() - started
            result = {
                "width_A": width,
                "wave_points": int(len(sme.wave[0])),
                "input_lines": int(len(linelist)),
                "parse_sec_excluded": parse_sec,
                "synthesis_wall_sec": wall,
                "flux_checksum": float(np.sum(np.asarray(sme.synth[0], dtype=float))),
            }
            output["models"].setdefault(name, {})[f"{width:g}"] = result
            args.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
            print(name, json.dumps(result, sort_keys=True), flush=True)
        del linelist
    return output


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--cache", type=Path, default=DEFAULT_CACHE)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument(
        "--widths", type=float, nargs="+", default=list(WIDTHS_A)
    )
    return parser.parse_args()


if __name__ == "__main__":
    ns = parse_args()
    results = run(ns)
    ns.output.write_text(json.dumps(results, indent=2, sort_keys=True) + "\n")
