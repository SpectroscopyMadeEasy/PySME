"""Validate native batched RKINTS against its prototype and legacy RKINTS."""

from __future__ import annotations

import argparse
import copy
import json
from pathlib import Path

import numpy as np

from almax_seeded_transfer_grid_benchmark import LINE_LIST, MODELS, make_sme
from batched_rkints_benchmark import (
    batched_refinement,
    legacy_initial_grid,
    normalized_intrinsic,
    raw_transf,
    restore_precomputed_line_state,
)
from pysme.linelist.vald import ValdFile
from pysme.synthesize import Synthesizer


ROOT = Path(__file__).resolve().parents[1]


def stats(delta):
    delta = np.asarray(delta, dtype=float)
    return {
        "max_abs": float(np.max(np.abs(delta))),
        "rms": float(np.sqrt(np.mean(delta * delta))),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "analysis" / "production_batched_rkints_validation.json",
    )
    args = parser.parse_args()
    linelist = ValdFile(str(LINE_LIST))
    result = {"models": {}}

    for name in MODELS:
        print(f"running {name}", flush=True)
        setup_wave = np.arange(5195.0, 5205.0 + 0.025, 0.05)
        sme = make_sme(name, copy.deepcopy(linelist), setup_wave)
        synth = Synthesizer()
        synth.synthesize_spectrum(
            sme,
            updateStructure=False,
            reuse_wavelength_grid=True,
            radial_velocity_mode="fast",
            linelist_mode="all",
            smelib_lineinfo_mode=2,
        )
        dll = synth.get_dll()
        restore_precomputed_line_state(dll, sme)
        legacy_wave, _, _, legacy_elapsed = raw_transf(
            dll, sme, None, keep_lineop=False, adaptive_mode="legacy"
        )

        supported = ~np.asarray(sme.line_ion_mask, dtype=bool)
        strong_all = np.asarray(sme.linelist["strong"], dtype=bool)
        line_wave = np.asarray(sme.linelist["wlcent"], dtype=float)
        initial_wave, _ = legacy_initial_grid(
            float(legacy_wave[0]),
            float(legacy_wave[-1]),
            line_wave,
            supported & strong_all,
        )
        restore_precomputed_line_state(dll, sme)
        ref_wave, ref_sint, ref_cint, _ = batched_refinement(
            dll, sme, initial_wave
        )
        ref_flux = normalized_intrinsic(sme, ref_sint, ref_cint)

        restore_precomputed_line_state(dll, sme)
        wave, sint, cint, elapsed = raw_transf(
            dll, sme, None, keep_lineop=False, adaptive_mode="batched"
        )
        flux = normalized_intrinsic(sme, sint, cint)
        expected_range = np.column_stack(
            (
                np.asarray(sme.linelist["line_range_s"], dtype=float)[supported],
                np.asarray(sme.linelist["line_range_e"], dtype=float)[supported],
            )
        )
        range_after = np.asarray(dll.GetLineRange(), dtype=float)
        result["models"][name] = {
            "points": int(wave.size),
            "seconds": elapsed,
            "legacy_points": int(legacy_wave.size),
            "legacy_seconds": legacy_elapsed,
            "speedup_vs_legacy": float(legacy_elapsed / elapsed),
            "physical_ranges_immutable": bool(
                np.array_equal(range_after, expected_range)
            ),
            "wavelength_exact": bool(np.array_equal(wave, ref_wave)),
            "wavelength_max_abs_A": float(np.max(np.abs(wave - ref_wave)))
            if wave.shape == ref_wave.shape
            else None,
            "flux_difference": stats(flux - ref_flux)
            if flux.shape == ref_flux.shape
            else None,
        }

    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(args.output)


if __name__ == "__main__":
    main()
