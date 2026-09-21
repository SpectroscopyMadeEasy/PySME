"""Benchmark the production adaptive SMElib continuum-opacity grid.

This script reuses the exact reference arrays produced by
``continuum_opacity_grid_benchmark.py`` and exercises the public production
API.  It writes a compact JSON report rather than modifying SMElib for
instrumentation.
"""

from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np

from pysme.linelist.linelist import LineList
from pysme.synthesize import Synthesizer

from continuum_opacity_grid_benchmark import (
    CDR_NPZ,
    FULL_LINELIST,
    MODELS,
    SMALL_LINELIST,
    make_sme,
    prepare_dll,
)


ROOT = Path(__file__).resolve().parents[1]
CONTINUUM_REFERENCE = ROOT / "analysis" / "continuum_opacity_grid_reference.npz"
OUTPUT = ROOT / "analysis" / "adaptive_continuum_grid_benchmark.json"
TOLERANCES = (3e-3, 1e-3, 3e-4, 1e-4)


def error_summary(approx, exact, floor=0.0):
    delta = np.abs(np.asarray(approx) - np.asarray(exact))
    scale = np.maximum(np.abs(exact), floor)
    relative = delta / np.maximum(scale, 1e-300)
    return {
        "max_abs": float(np.max(delta)),
        "rms_abs": float(np.sqrt(np.mean(delta * delta))),
        "p99_rel": float(np.percentile(relative, 99)),
        "max_rel": float(np.max(relative)),
    }


def spectrum_summary(approx, exact, wave):
    delta = np.asarray(approx) - np.asarray(exact)
    index = int(np.argmax(np.abs(delta)))
    return {
        "max_abs": float(np.max(np.abs(delta))),
        "rms": float(np.sqrt(np.mean(delta * delta))),
        "p99_abs": float(np.percentile(np.abs(delta), 99)),
        "worst_wavelength_A": float(wave[index]),
    }


def run_continuum_scan(output):
    with np.load(CONTINUUM_REFERENCE) as reference:
        wave = reference["wave"]
        exact = reference["chi"]

    sme = make_sme(*MODELS[0], SMALL_LINELIST)
    dll, _, _ = prepare_dll(sme, (3490.0, 9510.0))
    configs = {}
    for tolerance in TOLERANCES:
        dll.SetContinuumOpacityGrid(
            "adaptive", base_step=1.0, rtol=tolerance, min_step=1e-3
        )
        started = time.perf_counter()
        approximate = np.stack(
            [dll.GetContinuumOpacityComponents(float(item))[2] for item in wave]
        )
        elapsed = time.perf_counter() - started
        per_depth_floor = np.max(np.abs(exact), axis=1)[:, None] * 1e-12
        relative = np.abs(approximate - exact) / np.maximum(
            np.abs(exact), per_depth_floor
        )
        location = np.unravel_index(int(np.argmax(relative)), relative.shape)
        configs[f"{tolerance:g}"] = {
            "wall_sec": elapsed,
            "p99_relative": float(np.percentile(relative, 99)),
            "max_relative": float(relative[location]),
            "worst_wavelength_A": float(wave[location[0]]),
            "worst_depth_index": int(location[1]),
            "grid": dll.GetContinuumOpacityGridStats(),
        }
    output["continuum_scan"] = {
        "model": MODELS[0][0],
        "wavelength_count": int(wave.size),
        "wavelength_range_A": [float(wave.min()), float(wave.max())],
        "configs": configs,
    }


def install_line_selection(sme, central, range_s, range_e, strong):
    frame = sme.linelist._lines
    frame["central_depth"] = central
    frame["line_range_s"] = range_s
    frame["line_range_e"] = range_e
    frame["strong"] = strong
    sme.linelist.cdr_paras = np.array([sme.teff, sme.logg, sme.monh, sme.vmic])
    sme.linelist.cdr_paras_h_stark_convolution = "legacy"
    sme.linelist.cdr_paras_thres["strong_depth"] = 0.001
    sme.linelist.cdr_paras_thres["strong_bin_width"] = 0.2


def run_cdr(output):
    if not FULL_LINELIST.exists():
        raise FileNotFoundError(FULL_LINELIST)
    reference = np.load(CDR_NPZ)
    wave_spectrum = reference["wave"]
    models = {}
    calculated = {}
    for model in MODELS:
        name = model[0]
        print(f"adaptive CDR: {name}", flush=True)
        sme = make_sme(*model, FULL_LINELIST)
        dll, ion_mask, ionization_sec = prepare_dll(sme, (4998.0, 5011.0))
        keep = ~ion_mask
        wavelength = np.asarray(sme.linelist["wlcent"], dtype=float)
        species = np.char.strip(np.asarray(sme.linelist["species"], dtype=str))
        exact_central = reference[f"{name}_central"]
        exact_strong = reference[f"{name}_strong"].astype(bool)
        exact_range_s = reference[f"{name}_range_s"]
        exact_range_e = reference[f"{name}_range_e"]
        exact_spectrum = reference[f"{name}_spectrum"]
        model_configs = {}
        calculated[name] = {}

        for tolerance in TOLERANCES:
            dll.SetContinuumOpacityGrid(
                "adaptive", base_step=1.0, rtol=tolerance, min_step=1e-3
            )
            started = time.perf_counter()
            _, ranges_valid = dll.ALMAXRange(1e-4)
            almax_sec = time.perf_counter() - started
            started = time.perf_counter()
            central_valid = np.asarray(dll.CentralDepth(sme.mu, 1e-4), dtype=float)
            central_sec = time.perf_counter() - started

            central = np.full(wavelength.size, np.nan)
            central[keep] = central_valid
            central[species == "H 1"] = 1.0
            strong = np.asarray(
                Synthesizer.flag_strong_lines_by_bins(
                    wavelength,
                    central,
                    bin_width=0.2,
                    threshold=0.001,
                    valid_mask=keep,
                ),
                dtype=bool,
            )
            ranges_valid = np.asarray(ranges_valid, dtype=float)
            range_s = np.full(wavelength.size, np.nan)
            range_e = np.full(wavelength.size, np.nan)
            range_s[keep] = ranges_valid[:, 0]
            range_e[keep] = ranges_valid[:, 1]
            calculated[name][tolerance] = (central.copy(), strong.copy())

            finite = np.isfinite(central) & np.isfinite(exact_central)
            model_configs[f"{tolerance:g}"] = {
                "lineinfo_sec": almax_sec + central_sec,
                "almax_range_sec": almax_sec,
                "central_depth_sec": central_sec,
                "ionization_sec": ionization_sec,
                "grid": dll.GetContinuumOpacityGridStats(),
                "central_depth": error_summary(
                    central[finite], exact_central[finite], floor=1e-6
                ),
                "selection": {
                    "exact_selected": int(np.count_nonzero(exact_strong)),
                    "adaptive_selected": int(np.count_nonzero(strong)),
                    "false_positive": int(np.count_nonzero(strong & ~exact_strong)),
                    "false_negative": int(np.count_nonzero(exact_strong & ~strong)),
                },
                "range_finite_count": int(
                    np.count_nonzero(np.isfinite(range_s) & np.isfinite(range_e))
                ),
            }
        models[name] = {
            "atmosphere": dict(
                zip(("name", "teff", "logg", "monh", "vmic"), model)
            ),
            "configs": model_configs,
        }
        output["cdr"] = {"linelist": str(FULL_LINELIST), "models": models}
        OUTPUT.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")

    # Synthesizer may switch among per-thread native-library copies.  Keep the
    # direct low-level line-info measurements above in one phase, then perform
    # all high-level spectrum calls in a separate phase.
    for model in MODELS:
        name = model[0]
        print(f"adaptive spectra: {name}", flush=True)
        sme = make_sme(*model, FULL_LINELIST)
        exact_range_s = reference[f"{name}_range_s"]
        exact_range_e = reference[f"{name}_range_e"]
        exact_spectrum = reference[f"{name}_spectrum"]
        for tolerance in TOLERANCES:
            central, strong = calculated[name][tolerance]
            install_line_selection(
                sme, central, exact_range_s, exact_range_e, strong
            )
            sme.continuum_grid = "adaptive"
            sme.continuum_grid_base_step = 1.0
            sme.continuum_grid_rtol = tolerance
            sme.continuum_grid_min_step = 1e-3
            started = time.perf_counter()
            Synthesizer().synthesize_spectrum(
                sme, linelist_mode="dynamic", smelib_lineinfo_mode=2
            )
            synthesis_sec = time.perf_counter() - started
            spectrum = np.asarray(sme.synth[0], dtype=float).copy()
            config = models[name]["configs"][f"{tolerance:g}"]
            config["synthesis_sec"] = synthesis_sec
            config["spectrum"] = spectrum_summary(
                spectrum, exact_spectrum, wave_spectrum
            )
            OUTPUT.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
    reference.close()


def run_edge_spectrum(output):
    sme = make_sme(*MODELS[0], SMALL_LINELIST)
    linelist = LineList()
    linelist.add("Fe 1", 3756.62, 1.0, -2.0, 7.0, -6.0, 200.0)
    sme.linelist = linelist
    wave = np.linspace(3756.55, 3756.70, 301)
    sme.wave = [wave]
    sme.wint = [wave]
    sme.line_select_method = "internal"
    sme.line_select_policy = "auto"

    sme.continuum_grid = "exact"
    Synthesizer().synthesize_spectrum(sme)
    exact = np.asarray(sme.synth[0], dtype=float).copy()

    configs = {}
    for tolerance in TOLERANCES:
        sme.continuum_grid = "adaptive"
        sme.continuum_grid_rtol = tolerance
        Synthesizer().synthesize_spectrum(sme)
        approximate = np.asarray(sme.synth[0], dtype=float).copy()
        configs[f"{tolerance:g}"] = spectrum_summary(approximate, exact, wave)
    output["mg_edge_spectrum"] = {
        "model": MODELS[0][0],
        "range_A": [float(wave[0]), float(wave[-1])],
        "line_center_A": 3756.62,
        "configs": configs,
    }


def main():
    output = {
        "base_step_A": 1.0,
        "min_step_A": 1e-3,
        "tolerances": TOLERANCES,
        # One separately timed exact run of the identical ALMAXRange +
        # CentralDepth call chain. It is recorded here so normal benchmark
        # reruns do not pay the extra minute unless the baseline is refreshed.
        "exact_solar_same_path": {
            "ionization_sec": 0.19880191702395678,
            "almax_range_sec": 30.96389058290515,
            "central_depth_sec": 30.267129791900516,
            "lineinfo_sec": 61.23102037480567,
        },
    }
    run_continuum_scan(output)
    run_cdr(output)
    run_edge_spectrum(output)
    OUTPUT.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
    print(OUTPUT)


if __name__ == "__main__":
    main()
