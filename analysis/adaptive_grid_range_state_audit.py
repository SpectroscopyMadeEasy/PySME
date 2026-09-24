"""Isolate legacy RKINTS seed-scan range mutation from batched refinement."""

from __future__ import annotations

import argparse
import copy
import json
import os
from contextlib import contextmanager
from pathlib import Path

import numpy as np

from almax_seeded_transfer_grid_benchmark import LINE_LIST, MODELS, make_sme
from batched_rkints_benchmark import (
    grid_comparison,
    normalized_intrinsic,
    raw_transf,
    restore_precomputed_line_state,
)
from pysme.linelist.vald import ValdFile
from pysme.synthesize import Synthesizer


ROOT = Path(__file__).resolve().parents[1]


@contextmanager
def diagnostic_state(range_state):
    keys = ("SME_RKINTS_SECOND_PRUNING", "SME_RKINTS_RANGE_STATE")
    old = {key: os.environ.get(key) for key in keys}
    os.environ["SME_RKINTS_SECOND_PRUNING"] = "off"
    if range_state is None:
        os.environ.pop("SME_RKINTS_RANGE_STATE", None)
    else:
        os.environ["SME_RKINTS_RANGE_STATE"] = range_state
    try:
        yield
    finally:
        for key, value in old.items():
            if value is None:
                os.environ.pop(key, None)
            else:
                os.environ[key] = value


def error_stats(left, right):
    delta = np.asarray(left) - np.asarray(right)
    return {
        "max_abs": float(np.max(np.abs(delta))),
        "rms": float(np.sqrt(np.mean(delta * delta))),
        "mean": float(np.mean(delta)),
        "p99_abs": float(np.percentile(np.abs(delta), 99)),
    }


def compare(wave_a, flux_a, wave_b, flux_b):
    lo = max(float(wave_a[0]), float(wave_b[0]))
    hi = min(float(wave_a[-1]), float(wave_b[-1]))
    common = np.arange(lo, hi + 0.0003125, 0.000625)
    return {
        "grid": {
            key: value
            for key, value in grid_comparison(wave_a, wave_b).items()
            if not isinstance(value, np.ndarray)
        },
        "fine_resampled": error_stats(
            np.interp(common, wave_a, flux_a), np.interp(common, wave_b, flux_b)
        ),
    }


def run_model(name, linelist):
    setup_wave = np.arange(5195.0, 5205.0 + 0.025, 0.05)
    sme = make_sme(name, linelist, setup_wave)
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
    compact = ~np.asarray(sme.line_ion_mask, dtype=bool)
    baseline_ranges = np.column_stack(
        (
            np.asarray(sme.linelist["line_range_s"], dtype=float)[compact],
            np.asarray(sme.linelist["line_range_e"], dtype=float)[compact],
        )
    )

    restore_precomputed_line_state(dll, sme)
    batch_wave, batch_sint, batch_cint, batch_sec = raw_transf(
        dll, sme, None, keep_lineop=False, adaptive_mode="batched"
    )
    batch_flux = normalized_intrinsic(sme, batch_sint, batch_cint)

    restore_precomputed_line_state(dll, sme)
    with diagnostic_state(None):
        legacy_wave, legacy_sint, legacy_cint, legacy_sec = raw_transf(
            dll, sme, None, keep_lineop=False, adaptive_mode="legacy"
        )
    legacy_flux = normalized_intrinsic(sme, legacy_sint, legacy_cint)
    mutated_ranges = np.asarray(dll.GetLineRange(), dtype=float)
    range_changed = np.any(mutated_ranges != baseline_ranges, axis=1)

    restore_precomputed_line_state(dll, sme)
    with diagnostic_state("fixed"):
        fixed_wave, fixed_sint, fixed_cint, fixed_sec = raw_transf(
            dll, sme, None, keep_lineop=False, adaptive_mode="legacy"
        )
    fixed_flux = normalized_intrinsic(sme, fixed_sint, fixed_cint)

    changed_indices = np.flatnonzero(range_changed)
    samples = []
    for compact_index in changed_indices[:10]:
        python_index = int(np.flatnonzero(compact)[compact_index])
        samples.append(
            {
                "compact_index": int(compact_index),
                "python_index": python_index,
                "species": str(sme.linelist["species"][python_index]),
                "wavelength": float(sme.linelist["wlcent"][python_index]),
                "before": baseline_ranges[compact_index].tolist(),
                "after": mutated_ranges[compact_index].tolist(),
            }
        )

    return {
        "model": name,
        "stellar_parameters": dict(
            zip(("teff", "logg", "monh", "vmic"), MODELS[name])
        ),
        "points": {
            "batched": int(batch_wave.size),
            "legacy_mutable_ranges": int(legacy_wave.size),
            "legacy_immutable_ranges": int(fixed_wave.size),
        },
        "seconds": {
            "batched": batch_sec,
            "legacy_mutable_ranges": legacy_sec,
            "legacy_immutable_ranges": fixed_sec,
        },
        "range_mutation": {
            "changed_lines": int(changed_indices.size),
            "samples": samples,
        },
        "legacy_mutable_minus_batched": compare(
            legacy_wave, legacy_flux, batch_wave, batch_flux
        ),
        "legacy_immutable_minus_batched": compare(
            fixed_wave, fixed_flux, batch_wave, batch_flux
        ),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--models", nargs="+", default=["solar_dwarf", "giant"], choices=MODELS
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "analysis" / "adaptive_grid_range_state_audit.json",
    )
    args = parser.parse_args()
    linelist = ValdFile(str(LINE_LIST))
    result = {
        "setup": {
            "line_list": str(LINE_LIST),
            "window_A": [5195.0, 5205.0],
            "second_pruning": "disabled in both legacy diagnostic runs",
            "diagnostic_switch": "SME_RKINTS_RANGE_STATE=fixed",
        },
        "models": {},
    }
    for name in args.models:
        print(f"running {name}", flush=True)
        result["models"][name] = run_model(name, copy.deepcopy(linelist))
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(args.output)


if __name__ == "__main__":
    main()
