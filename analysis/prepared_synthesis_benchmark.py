"""Benchmark repeated abundance synthesis with persistent SMElib state.

This is a bounded prototype for a possible ``PreparedSynthesis`` lifecycle.
It compares the current repeated-synthesis path, which re-inputs the line list
and model on every call, with a prepared path that keeps those unchanged inputs
resident in SMElib.  Both paths still run InputAbund, Ionization, Opacity, and
Transf for every abundance point and both use the exact EOS warm-start mode.

The initial atmosphere interpolation and ALMAX/range preparation are excluded
from timed evaluations in both paths.  No public API is changed by this script.
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

from line_centric_wide_benchmark import (
    DEFAULT_CACHE,
    DEFAULT_SOURCE,
    MODELS,
    STEP_A,
    ensure_window,
    _window_bounds,
)


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = ROOT / "analysis" / "prepared_synthesis_benchmark.json"
DEFAULT_WIDTHS_A = (10.0, 200.0)
FE_OFFSETS_DEX = (0.00, 0.05, -0.05, 0.10, -0.10, 0.00)


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
    sme.line_select_recompute = "if_stale"
    sme.continuum_grid = "adaptive"
    sme.continuum_grid_base_step = 1.0
    sme.continuum_grid_rtol = 1e-3
    sme.continuum_grid_min_step = 1e-3
    return sme


def synthesize(synth: Synthesizer, sme: SME_Structure, **kwargs):
    started = time.perf_counter()
    result = synth.synthesize_spectrum(
        sme,
        updateStructure=False,
        reuse_wavelength_grid=True,
        radial_velocity_mode="fast",
        linelist_mode="all",
        smelib_lineinfo_mode=0,
        **kwargs,
    )
    wall = time.perf_counter() - started
    return wall, np.asarray(result[1][0], dtype=np.float64).copy()


def reset_eos_warm_start(dll) -> None:
    # Each setter call clears native EOS history, making the two A/B arms start
    # from the same empty-history state before their untimed preparation call.
    dll.SetEosWarmStartMode(False)
    dll.SetEosWarmStartMode(True)


def run_arm(
    mode: str,
    model: tuple,
    linelist: ValdFile,
    width: float,
) -> tuple[dict, list[np.ndarray]]:
    sme = make_sme(model, linelist, width)
    synth = Synthesizer()
    dll = synth.get_dll()
    reset_eos_warm_start(dll)

    base_fe = float(sme.abund.A["Fe"])
    sme.abund.A["Fe"] = base_fe
    prepare_wall, _ = synthesize(
        synth,
        sme,
        passLineList=True,
        passAtmosphere=True,
    )

    times = []
    spectra = []
    for offset in FE_OFFSETS_DEX:
        sme.abund.A["Fe"] = base_fe + offset
        if mode == "current":
            wall, flux = synthesize(
                synth,
                sme,
                passLineList=True,
                passAtmosphere=True,
            )
        elif mode == "prepared":
            wall, flux = synthesize(
                synth,
                sme,
                passLineList=False,
                passAtmosphere=False,
                passAbund=True,
            )
        else:
            raise ValueError(mode)
        times.append(wall)
        spectra.append(flux)
        print(
            f"{model[0]} width={width:g} mode={mode} "
            f"Fe_offset={offset:+.2f} wall={wall:.6f}s",
            flush=True,
        )

    dll.SetEosWarmStartMode(False)
    values = np.asarray(times, dtype=float)
    record = {
        "untimed_prepare_sec": prepare_wall,
        "evaluation_sec": values.tolist(),
        "total_evaluation_sec": float(np.sum(values)),
        "median_evaluation_sec": float(np.median(values)),
        "min_evaluation_sec": float(np.min(values)),
    }
    return record, spectra


def run(args: argparse.Namespace) -> dict:
    output = {
        "description": (
            "Repeated Fe-abundance synthesis: current full re-input versus "
            "persistent line-list/model state"
        ),
        "source": str(args.source),
        "step_A": STEP_A,
        "fe_offsets_dex": list(FE_OFFSETS_DEX),
        "models": {},
    }
    for width in args.widths:
        path = ensure_window(args.source, args.cache, width)
        started = time.perf_counter()
        linelist = ValdFile(str(path))
        print(
            f"loaded width={width:g}: {len(linelist):,} lines in "
            f"{time.perf_counter() - started:.2f}s",
            flush=True,
        )
        for model in MODELS:
            current, current_fluxes = run_arm("current", model, linelist, width)
            prepared, prepared_fluxes = run_arm("prepared", model, linelist, width)
            max_abs = []
            rms = []
            for reference, candidate in zip(current_fluxes, prepared_fluxes):
                delta = candidate - reference
                max_abs.append(float(np.max(np.abs(delta))))
                rms.append(float(np.sqrt(np.mean(delta * delta))))
            comparison = {
                "width_A": width,
                "wave_points": int(len(current_fluxes[0])),
                "input_lines": int(len(linelist)),
                "current": current,
                "prepared": prepared,
                "total_speedup": (
                    current["total_evaluation_sec"]
                    / prepared["total_evaluation_sec"]
                ),
                "median_speedup": (
                    current["median_evaluation_sec"]
                    / prepared["median_evaluation_sec"]
                ),
                "max_abs_flux_delta_by_evaluation": max_abs,
                "rms_flux_delta_by_evaluation": rms,
                "max_abs_flux_delta": float(max(max_abs)),
                "max_rms_flux_delta": float(max(rms)),
            }
            output["models"].setdefault(model[0], {})[f"{width:g}"] = comparison
            args.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
            print(json.dumps(comparison, sort_keys=True), flush=True)
        del linelist
    return output


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--cache", type=Path, default=DEFAULT_CACHE)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument(
        "--widths", type=float, nargs="+", default=list(DEFAULT_WIDTHS_A)
    )
    return parser.parse_args()


if __name__ == "__main__":
    ns = parse_args()
    results = run(ns)
    ns.output.write_text(json.dumps(results, indent=2, sort_keys=True) + "\n")
