"""Summarize deterministic spectrum/line-state differences in the v1.2 audit."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "analysis" / "v120_release_benchmark.json"
OUTPUT = ROOT / "analysis" / "v120_release_accuracy.json"


def load_first(profile: dict) -> dict[str, np.ndarray]:
    path = ROOT / profile["runs"][0]["array_file"]
    with np.load(path) as data:
        return {name: np.asarray(data[name]) for name in data.files}


def compare(reference: dict[str, np.ndarray], candidate: dict[str, np.ndarray]) -> dict:
    if not np.array_equal(reference["wave"], candidate["wave"]):
        raise ValueError("Post-integration wavelength grids differ")
    wave = reference["wave"]
    delta = candidate["flux"] - reference["flux"]
    absolute = np.abs(delta)
    reference_ew = float(np.trapezoid(1.0 - reference["flux"], wave))
    candidate_ew = float(np.trapezoid(1.0 - candidate["flux"], wave))
    result = {
        "bitwise_flux_equal": bool(np.array_equal(reference["flux"], candidate["flux"])),
        "max_abs_flux": float(np.max(absolute)),
        "rms_flux": float(np.sqrt(np.mean(delta * delta))),
        "p99_abs_flux": float(np.percentile(absolute, 99)),
        "reference_ew_A": reference_ew,
        "candidate_ew_A": candidate_ew,
        "delta_ew_A": candidate_ew - reference_ew,
        "relative_ew": (
            (candidate_ew - reference_ew) / reference_ew
            if reference_ew != 0.0
            else None
        ),
    }

    if np.any(np.isfinite(reference["almax"])) and np.any(
        np.isfinite(candidate["almax"])
    ):
        almax_delta = candidate["almax"] - reference["almax"]
        finite = np.isfinite(almax_delta)
        result["almax"] = {
            "bitwise_equal": bool(np.array_equal(reference["almax"], candidate["almax"])),
            "max_abs": float(np.max(np.abs(almax_delta[finite]))),
        }
    if reference["strong"].shape == candidate["strong"].shape:
        removed = reference["strong"] & ~candidate["strong"]
        added = candidate["strong"] & ~reference["strong"]
        result["line_selection"] = {
            "reference_count": int(np.count_nonzero(reference["strong"])),
            "candidate_count": int(np.count_nonzero(candidate["strong"])),
            "removed": int(np.count_nonzero(removed)),
            "added": int(np.count_nonzero(added)),
        }
    if np.any(np.isfinite(reference["range_s"])) and np.any(
        np.isfinite(candidate["range_s"])
    ):
        range_delta = np.concatenate(
            [
                candidate["range_s"] - reference["range_s"],
                candidate["range_e"] - reference["range_e"],
            ]
        )
        finite = np.isfinite(range_delta)
        result["physical_ranges"] = {
            "bitwise_equal": bool(
                np.array_equal(reference["range_s"], candidate["range_s"])
                and np.array_equal(reference["range_e"], candidate["range_e"])
            ),
            "max_abs_A": float(np.max(np.abs(range_delta[finite]))),
            "changed_endpoints": int(np.count_nonzero(range_delta[finite])),
        }
    return result


def main():
    benchmark = json.loads(INPUT.read_text())
    output = {"comparisons": {}}
    pairs = {
        "continuum_only_C0B0_vs_C1B0": ("C0B0", "C1B0"),
        "batching_only_C0B0_vs_C0B1": ("C0B0", "C0B1"),
        "combined_C0B0_vs_C1B1": ("C0B0", "C1B1"),
        "historical_H_vs_C1B1": ("H", "C1B1"),
    }
    for case_name, case in benchmark["cases"].items():
        profiles = case["profiles"]
        comparisons = {}
        for name, (reference_name, candidate_name) in pairs.items():
            if reference_name not in profiles or candidate_name not in profiles:
                continue
            comparisons[name] = compare(
                load_first(profiles[reference_name]),
                load_first(profiles[candidate_name]),
            )
        if comparisons:
            output["comparisons"][case_name] = comparisons
    OUTPUT.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
    print(OUTPUT)


if __name__ == "__main__":
    main()
