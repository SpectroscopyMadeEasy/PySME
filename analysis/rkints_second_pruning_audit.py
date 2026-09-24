"""Audit RKINTS' historical disk-centre second weak-line pruning.

This is a diagnostic benchmark.  It relies on the temporary SMElib environment
controls ``SME_RKINTS_SECOND_PRUNING`` and ``SME_RKINTS_PRUNE_LOG`` and never
changes the production default.  The main comparison uses the same stellar
models and 5195--5205 A window as ``batched_rkints_benchmark.py``.
"""

from __future__ import annotations

import argparse
import contextlib
import copy
import csv
import json
import os
import re
import tempfile
import time
from pathlib import Path

import numpy as np
from scipy.ndimage import gaussian_filter1d

from almax_seeded_transfer_grid_benchmark import (
    HI,
    LINE_LIST,
    LO,
    MODELS,
    equivalent_widths,
    error_stats,
    feature_windows,
    integrate_mu_areas,
    is_molecular,
    make_sme,
)
from batched_rkints_benchmark import raw_transf
from pysme.linelist.vald import ValdFile
from pysme.synthesize import Synthesizer


ROOT = Path(__file__).resolve().parents[1]
BASE_BATCH_NPZ = ROOT / "analysis" / "batched_rkints_benchmark.npz"
FINE_STEP = 0.000625
LOW_RANGE_FLOOR = 1e-6
# SMElib reads this switch once on the first Transf call.  Enable it before
# benchmark setup so later per-call file-descriptor captures contain native
# OPMTRX/TBINTG/interval counters.
os.environ.setdefault("SME_TIMING", "1")
ELEMENTS = {
    "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg",
    "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr",
    "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br",
    "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd",
    "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La",
    "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er",
    "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au",
    "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th",
    "Pa", "U",
}
FE_GROUP = {"Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni"}


def json_clean(value):
    if isinstance(value, dict):
        return {str(key): json_clean(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_clean(item) for item in value]
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    return value


def normalized_intrinsic(sme, sint, cint):
    line = integrate_mu_areas(sme.mu, sint)
    continuum = integrate_mu_areas(sme.mu, cint)
    return np.asarray(line / continuum, dtype=float)


def compact_indices(sme):
    return np.flatnonzero(~np.asarray(sme.line_ion_mask, dtype=bool))


def input_line_state(dll, sme, ranges, mask):
    kept = compact_indices(sme)
    ranges = np.asarray(ranges, dtype=float)
    mask = np.asarray(mask, dtype=bool)
    if ranges.shape == (len(sme.linelist), 2):
        ranges = ranges[kept]
    if mask.shape == (len(sme.linelist),):
        mask = mask[kept]
    if ranges.shape != (kept.size, 2) or mask.shape != (kept.size,):
        raise ValueError("line-info arrays do not match compact SMElib line list")
    dll.SetLineInfoMode(2)
    dll.InputLinePrecomputedInfo(ranges[:, 0], ranges[:, 1], mask.astype(np.uint8))


@contextlib.contextmanager
def capture_native_stderr(path):
    """Capture C/C++ writes to file descriptor 2, not only Python stderr."""
    saved = os.dup(2)
    target = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_TRUNC, 0o600)
    try:
        os.dup2(target, 2)
        yield
    finally:
        os.dup2(saved, 2)
        os.close(saved)
        os.close(target)


def parse_native_timing(text):
    blocks = text.split("SMElib timing [Transf#")
    if len(blocks) < 2:
        return {}
    block = blocks[-1]

    def number(pattern, cast=float, default=None):
        match = re.search(pattern, block)
        return cast(match.group(1)) if match else default

    return {
        "total_sec": number(r"total\s+:\s+([0-9.]+)"),
        "setup_sec": number(r"setup\s+:\s+([0-9.]+)"),
        "lineopac_sec": number(r"LINEOPAC\s+:\s+([0-9.]+)"),
        "lineopac_calls": number(r"LINEOPAC\s+:.*calls=(\d+)", int),
        "opmtrx_sec": number(r"OPMTRX\s+:\s+([0-9.]+)"),
        "opmtrx_calls": number(r"OPMTRX\s+:.*calls=(\d+)", int),
        "tbintg_sec": number(r"TBINTG\s+:\s+([0-9.]+)"),
        "tbintg_calls": number(r"TBINTG\s+:.*calls=(\d+)", int),
        "interval_calls": number(r"intervals\s+: calls=(\d+)", int),
        "interval_candidates": number(r"intervals\s+: calls=\d+, candidates=(\d+)", int),
        "interval_full_scan": number(r"intervals\s+: calls=\d+, candidates=\d+, full_scan=(\d+)", int),
        "interval_reduction_percent": number(r"reduction=([0-9.]+)%"),
        "active_lines": number(r"active\(mark=0\)=(\d+)", int),
        "inactive_lines": number(r"inactive\(mark!=0\)=(\d+)", int),
    }


def timed_fixed_transf(dll, sme, wave, ranges, mask):
    input_line_state(dll, sme, ranges, mask)
    previous = os.environ.get("SME_TIMING")
    os.environ["SME_TIMING"] = "1"
    with tempfile.NamedTemporaryFile(prefix="rkints-timing-", suffix=".log", delete=False) as handle:
        timing_path = handle.name
    try:
        with capture_native_stderr(timing_path):
            got_wave, sint, cint, wall = raw_transf(
                dll, sme, np.asarray(wave, dtype=float), keep_lineop=False
            )
        text = Path(timing_path).read_text()
    finally:
        Path(timing_path).unlink(missing_ok=True)
        if previous is None:
            os.environ.pop("SME_TIMING", None)
        else:
            os.environ["SME_TIMING"] = previous
    native = parse_native_timing(text)
    native["python_wall_sec"] = wall
    return got_wave, normalized_intrinsic(sme, sint, cint), native


def median_timing(records):
    keys = sorted({key for record in records for key in record})
    output = {"repeats": len(records)}
    for key in keys:
        values = [record[key] for record in records if record.get(key) is not None]
        if values:
            output[key] = float(np.median(values))
    return output


def read_prune_log(path, sme):
    kept = compact_indices(sme)
    rows = []
    with open(path, newline="") as handle:
        for raw in csv.DictReader(handle):
            compact = int(raw["line"])
            full = int(kept[compact])
            rows.append(
                {
                    "phase": raw["phase"],
                    "sequence": int(raw["sequence"]),
                    "compact_index": compact,
                    "line_index": full,
                    "wavelength": float(raw["wavelength"]),
                    "disk_intensity": float(raw["disk_intensity"]),
                    "continuum_intensity": float(raw["continuum_intensity"]),
                    "depression": float(raw["depression"]),
                    "would_prune": bool(int(raw["would_prune"])),
                    "applied": bool(int(raw["applied"])),
                    "mode": int(raw["mode"]),
                }
            )
    return rows


def run_prune_diagnostic(dll, sme, mode, log_path):
    prod_ranges = np.column_stack(
        (
            np.asarray(sme.linelist["line_range_s"], dtype=float),
            np.asarray(sme.linelist["line_range_e"], dtype=float),
        )
    )
    prod_mask = np.asarray(sme.linelist["strong"], dtype=bool)
    input_line_state(dll, sme, prod_ranges, prod_mask)
    old_mode = os.environ.get("SME_RKINTS_SECOND_PRUNING")
    old_log = os.environ.get("SME_RKINTS_PRUNE_LOG")
    os.environ["SME_RKINTS_SECOND_PRUNING"] = mode
    os.environ["SME_RKINTS_PRUNE_LOG"] = str(log_path)
    try:
        wave, sint, cint, elapsed = raw_transf(
            dll, sme, None, keep_lineop=False, adaptive_mode="legacy"
        )
    finally:
        if old_mode is None:
            os.environ.pop("SME_RKINTS_SECOND_PRUNING", None)
        else:
            os.environ["SME_RKINTS_SECOND_PRUNING"] = old_mode
        if old_log is None:
            os.environ.pop("SME_RKINTS_PRUNE_LOG", None)
        else:
            os.environ["SME_RKINTS_PRUNE_LOG"] = old_log
    rows = read_prune_log(log_path, sme)
    return wave, normalized_intrinsic(sme, sint, cint), elapsed, rows


def species_category(label):
    label = str(label).strip()
    token = label.split()[0] if label else ""
    if token.startswith("TiO"):
        return "TiO"
    if token.startswith("MgH"):
        return "MgH"
    if token in FE_GROUP:
        return "Fe-group"
    if token in ELEMENTS:
        return "other atomic"
    return "other molecular"


def almax_bin(value):
    if value < 1e-3:
        return "1e-4--1e-3"
    if value < 1e-2:
        return "1e-3--1e-2"
    if value < 1e-1:
        return "1e-2--1e-1"
    return ">=1e-1"


def enrich_removed_rows(rows, sme):
    almax = np.asarray(sme.linelist["almax_ratio"], dtype=float)
    range_s = np.asarray(sme.linelist["line_range_s"], dtype=float)
    range_e = np.asarray(sme.linelist["line_range_e"], dtype=float)
    species = np.asarray(sme.linelist["species"]).astype(str)
    output = []
    for row in rows:
        if row["phase"] != "seed" or not row["applied"]:
            continue
        index = row["line_index"]
        item = dict(row)
        item.update(
            {
                "species": species[index],
                "category": species_category(species[index]),
                "almax": float(almax[index]),
                "almax_bin": almax_bin(float(almax[index])),
                "range_start": float(range_s[index]),
                "range_end": float(range_e[index]),
            }
        )
        output.append(item)
    return output


def count_by(items, key):
    result = {}
    for item in items:
        label = str(item[key])
        result[label] = result.get(label, 0) + 1
    return dict(sorted(result.items()))


def direction_stats(delta):
    delta = np.asarray(delta, dtype=float)
    significant = np.abs(delta) > 1e-14
    if np.any(significant):
        positive_nonzero = float(np.mean(delta[significant] > 0))
        negative_nonzero = float(np.mean(delta[significant] < 0))
    else:
        positive_nonzero = negative_nonzero = 0.0
    return {
        "fraction_positive_all": float(np.mean(delta > 0)),
        "fraction_negative_all": float(np.mean(delta < 0)),
        "fraction_positive_nonzero": positive_nonzero,
        "fraction_negative_nonzero": negative_nonzero,
        "mean": float(np.mean(delta)),
        "median": float(np.median(delta)),
        "minimum": float(np.min(delta)),
        "maximum": float(np.max(delta)),
    }


def metrics_on_common_grid(candidate, reference, wave, windows):
    result = {
        "intrinsic": error_stats(candidate, reference),
        "post_convolution": {},
        "equivalent_width_delta_A": {},
        "signed_difference": direction_stats(candidate - reference),
    }
    step = float(wave[1] - wave[0])
    for resolution in (20000.0, 60000.0):
        sigma_pixels = (float(np.mean(wave)) / resolution) / 2.354820045 / step
        result["post_convolution"][str(int(resolution))] = error_stats(
            gaussian_filter1d(candidate, sigma_pixels, mode="nearest"),
            gaussian_filter1d(reference, sigma_pixels, mode="nearest"),
        )
    candidate_ew = equivalent_widths(wave, candidate, windows)
    reference_ew = equivalent_widths(wave, reference, windows)
    result["equivalent_width_delta_A"] = {
        name: candidate_ew[name] - reference_ew[name] for name in windows
    }
    return result


def setup_model(model_name, linelist):
    sme = make_sme(model_name, linelist, np.arange(LO, HI + 0.025, 0.05))
    synth = Synthesizer()
    started = time.perf_counter()
    synth.synthesize_spectrum(
        sme,
        updateStructure=False,
        reuse_wavelength_grid=True,
        radial_velocity_mode="fast",
        linelist_mode="all",
        smelib_lineinfo_mode=2,
    )
    return sme, synth, time.perf_counter() - started


def setup_window_model(model_name, linelist, lo, hi):
    sme = make_sme(model_name, linelist, np.arange(lo, hi + 0.025, 0.05))
    output_wave = np.arange(lo, hi + 0.025, 0.05)
    sme.wave = [output_wave]
    sme.wint = [output_wave]
    synth = Synthesizer()
    started = time.perf_counter()
    synth.synthesize_spectrum(
        sme,
        updateStructure=False,
        reuse_wavelength_grid=True,
        radial_velocity_mode="fast",
        linelist_mode="all",
        smelib_lineinfo_mode=2,
    )
    return sme, synth, time.perf_counter() - started


def compute_low_floor_lineinfo(dll, sme, wfirst, wlast, floor=LOW_RANGE_FLOOR):
    dll.InputWaveRange(float(wfirst), float(wlast))
    dll.Opacity()
    almax_compact, ranges_compact = dll.ALMAXRange(accrt=floor)
    kept = compact_indices(sme)
    almax = np.full(len(sme.linelist), np.nan, dtype=float)
    ranges = np.column_stack(
        (
            np.asarray(sme.linelist["wlcent"], dtype=float),
            np.asarray(sme.linelist["wlcent"], dtype=float),
        )
    )
    almax[kept] = np.asarray(almax_compact, dtype=float)
    ranges[kept] = np.asarray(ranges_compact, dtype=float)
    return almax, ranges


def run_performance_pair(dll, sme, grid, ranges, mask_l, mask_a, repeats=3):
    records = {"L_final_mask": [], "A_almax_only": []}
    last_flux = {}
    for _ in range(repeats):
        for name, mask in (("A_almax_only", mask_a), ("L_final_mask", mask_l)):
            _, flux, timing = timed_fixed_transf(dll, sme, grid, ranges, mask)
            records[name].append(timing)
            last_flux[name] = flux
    return {name: median_timing(value) for name, value in records.items()}, last_flux


def run_order_audit(dll, sme, workdir, fine_wave, prod_ranges, prod_mask):
    result = {}
    masks = {}
    fluxes = {}
    for mode in ("audit-forward", "audit-reverse"):
        log = workdir / f"order-{mode}.csv"
        _, _, elapsed, rows = run_prune_diagnostic(dll, sme, mode, log)
        decisions = [row for row in rows if row["phase"] == "decision" and row["applied"]]
        removed = np.asarray([row["line_index"] for row in decisions], dtype=int)
        mask = prod_mask.copy()
        mask[removed] = False
        masks[mode] = mask
        _, flux, _ = timed_fixed_transf(dll, sme, fine_wave, prod_ranges, mask)
        fluxes[mode] = flux
        result[mode] = {
            "diagnostic_transf_sec": elapsed,
            "removed_count": int(removed.size),
            "removed_indices": removed,
        }
    forward = set(result["audit-forward"]["removed_indices"].tolist())
    reverse = set(result["audit-reverse"]["removed_indices"].tolist())
    result["comparison"] = {
        "intersection": int(len(forward & reverse)),
        "forward_only": int(len(forward - reverse)),
        "reverse_only": int(len(reverse - forward)),
        "symmetric_difference": int(len(forward ^ reverse)),
        "spectrum_reverse_minus_forward": error_stats(
            fluxes["audit-reverse"], fluxes["audit-forward"]
        ),
    }
    return result, fluxes


def run_remove_one_audit(
    dll,
    sme,
    fine_wave,
    prod_ranges,
    prod_mask,
    removed_rows,
    flux_a,
    flux_l_final,
    maximum_difference_wavelength,
):
    candidates = [
        row
        for row in removed_rows
        if row["range_start"] < maximum_difference_wavelength < row["range_end"]
        or abs(row["wavelength"] - maximum_difference_wavelength) < 0.5
    ]
    candidates = sorted(
        candidates,
        key=lambda row: (-row["almax"], abs(row["wavelength"] - maximum_difference_wavelength)),
    )[:24]
    marginals = []
    for row in candidates:
        mask = prod_mask.copy()
        mask[row["line_index"]] = False
        _, flux, _ = timed_fixed_transf(dll, sme, fine_wave, prod_ranges, mask)
        marginal = np.abs(flux - flux_a)
        item = dict(row)
        item["marginal_max_abs"] = float(np.max(marginal))
        item["marginal_at_target"] = float(
            flux[np.argmin(np.abs(fine_wave - maximum_difference_wavelength))]
            - flux_a[np.argmin(np.abs(fine_wave - maximum_difference_wavelength))]
        )
        marginals.append(item)
    marginals.sort(key=lambda row: row["marginal_max_abs"], reverse=True)

    cumulative = {}
    for count in (1, 5, 20):
        chosen = marginals[: min(count, len(marginals))]
        mask = prod_mask.copy()
        mask[[row["line_index"] for row in chosen]] = False
        _, flux, _ = timed_fixed_transf(dll, sme, fine_wave, prod_ranges, mask)
        cumulative[str(count)] = {
            "actual_count": int(len(chosen)),
            "toward_final_removed_mask": error_stats(flux, flux_l_final),
            "change_from_almax_only": error_stats(flux, flux_a),
            "line_indices": [row["line_index"] for row in chosen],
        }
    cumulative["all"] = {
        "actual_count": int(len(removed_rows)),
        "toward_final_removed_mask": error_stats(flux_l_final, flux_l_final),
        "change_from_almax_only": error_stats(flux_l_final, flux_a),
        "line_indices": [row["line_index"] for row in removed_rows],
    }
    return {"candidate_marginals": marginals, "cumulative": cumulative}


def run_model(model_name, linelist, workdir, base_arrays, do_order=False, do_remove=False):
    sme, synth, setup_sec = setup_model(model_name, copy.deepcopy(linelist))
    dll = synth.get_dll()
    prod_mask = np.asarray(sme.linelist["strong"], dtype=bool)
    prod_ranges = np.column_stack(
        (
            np.asarray(sme.linelist["line_range_s"], dtype=float),
            np.asarray(sme.linelist["line_range_e"], dtype=float),
        )
    )
    supported = ~np.asarray(sme.line_ion_mask, dtype=bool)

    prune_log = workdir / f"{model_name}-immediate.csv"
    legacy_wave, legacy_flux, legacy_sec, prune_rows = run_prune_diagnostic(
        dll, sme, "immediate", prune_log
    )
    removed_rows = enrich_removed_rows(prune_rows, sme)
    removed_indices = np.asarray([row["line_index"] for row in removed_rows], dtype=int)
    legacy_final_mask = prod_mask.copy()
    legacy_final_mask[removed_indices] = False

    batch_grid = np.asarray(base_arrays[f"{model_name}_batch_wave"], dtype=float)
    performance, perf_flux = run_performance_pair(
        dll, sme, batch_grid, prod_ranges, legacy_final_mask, prod_mask
    )

    wfirst = float(legacy_wave[0])
    wlast = float(legacy_wave[-1])
    fine_wave = np.arange(wfirst, wlast + 0.5 * FINE_STEP, FINE_STEP)
    low_almax, low_ranges = compute_low_floor_lineinfo(dll, sme, wfirst, wlast)
    valid_low = supported & np.isfinite(low_almax)
    masks = {
        "L_final_mask": legacy_final_mask,
        "A_almax_only": prod_mask,
        "P_1e-5": valid_low & (low_almax >= 1e-5),
        "F_1e-6": valid_low & (low_almax >= LOW_RANGE_FLOOR),
    }
    ranges = {
        "L_final_mask": prod_ranges,
        "A_almax_only": prod_ranges,
        "P_1e-5": low_ranges,
        "F_1e-6": low_ranges,
    }
    fine_fluxes = {}
    fine_timings = {}
    for name in masks:
        _, flux, timing = timed_fixed_transf(
            dll, sme, fine_wave, ranges[name], masks[name]
        )
        fine_fluxes[name] = flux
        fine_timings[name] = timing

    centre = (fine_wave >= LO) & (fine_wave <= HI)
    compare_wave = fine_wave[centre]
    reference = fine_fluxes["F_1e-6"][centre]
    line_wave = np.asarray(sme.linelist["wlcent"], dtype=float)
    almax_prod = np.asarray(sme.linelist["almax_ratio"], dtype=float)
    molecule = is_molecular(np.asarray(sme.linelist["species"]))
    windows = feature_windows(line_wave, almax_prod, molecule)

    # Historical L and generation-batched A retain their own adaptive grids.
    batch_wave = np.asarray(base_arrays[f"{model_name}_batch_wave"], dtype=float)
    batch_flux = np.asarray(base_arrays[f"{model_name}_batch_flux"], dtype=float)
    adaptive = {
        "L_historical": np.interp(compare_wave, legacy_wave, legacy_flux),
        "A_batched": np.interp(compare_wave, batch_wave, batch_flux),
    }
    comparisons = {
        name: metrics_on_common_grid(flux, reference, compare_wave, windows)
        for name, flux in adaptive.items()
    }
    for name, flux in fine_fluxes.items():
        comparisons[name] = metrics_on_common_grid(
            flux[centre], reference, compare_wave, windows
        )

    direct_la = metrics_on_common_grid(
        adaptive["L_historical"], adaptive["A_batched"], compare_wave, windows
    )
    direct_la["definition"] = "L_historical minus A_batched"
    direct_la["signed_difference"] = direction_stats(
        adaptive["L_historical"] - adaptive["A_batched"]
    )

    category_counts = count_by(removed_rows, "category")
    species_counts = count_by(removed_rows, "species")
    bin_counts = count_by(removed_rows, "almax_bin")
    seed_rows = [row for row in prune_rows if row["phase"] == "seed"]

    result = {
        "model": model_name,
        "stellar_parameters": dict(zip(("teff", "logg", "monh", "vmic"), MODELS[model_name])),
        "setup_sec": setup_sec,
        "line_counts": {
            "input_legal": int(np.sum(supported)),
            "after_almax": int(np.sum(prod_mask)),
            "almax_overlapping_transfer_window": int(
                np.sum(
                    prod_mask
                    & (prod_ranges[:, 1] > wfirst)
                    & (prod_ranges[:, 0] < wlast)
                )
            ),
            "seed_centres_tested": int(len(seed_rows)),
            "additionally_mark2": int(len(removed_rows)),
            "fraction_of_almax_active": float(len(removed_rows) / max(1, np.sum(prod_mask))),
            "fraction_of_overlapping_almax": float(
                len(removed_rows)
                / max(
                    1,
                    np.sum(
                        prod_mask
                        & (prod_ranges[:, 1] > wfirst)
                        & (prod_ranges[:, 0] < wlast)
                    ),
                )
            ),
            "fraction_of_tested_centres": float(len(removed_rows) / max(1, len(seed_rows))),
            "permissive_1e-5": int(np.sum(masks["P_1e-5"])),
            "full_1e-6": int(np.sum(masks["F_1e-6"])),
        },
        "removed_summary": {
            "by_category": category_counts,
            "by_species": species_counts,
            "by_almax_bin": bin_counts,
            "maximum_removed_almax": float(max((row["almax"] for row in removed_rows), default=0.0)),
            "rows": removed_rows,
        },
        "legacy_adaptive": {
            "points": int(legacy_wave.size),
            "transf_sec": legacy_sec,
        },
        "same_grid_performance": performance,
        "fine_grid": {
            "wavelength_start": wfirst,
            "wavelength_end": wlast,
            "step_A": FINE_STEP,
            "points": int(fine_wave.size),
            "mode_line_counts": {name: int(np.sum(mask)) for name, mask in masks.items()},
            "timings": fine_timings,
        },
        "comparisons_to_F_1e-6": comparisons,
        "historical_L_minus_batched_A": direct_la,
        "feature_windows": {name: list(bounds) for name, bounds in windows.items()},
    }
    arrays = {
        f"{model_name}_legacy_wave": legacy_wave,
        f"{model_name}_legacy_flux": legacy_flux,
        f"{model_name}_fine_wave": fine_wave,
        f"{model_name}_low_almax": low_almax,
    }
    for name, flux in fine_fluxes.items():
        arrays[f"{model_name}_{name}_fine_flux"] = flux
    arrays[f"{model_name}_L_historical_compare_flux"] = adaptive["L_historical"]
    arrays[f"{model_name}_A_batched_compare_flux"] = adaptive["A_batched"]

    if do_order:
        order_result, order_fluxes = run_order_audit(
            dll, sme, workdir, fine_wave, prod_ranges, prod_mask
        )
        result["order_audit"] = order_result
        for name, flux in order_fluxes.items():
            arrays[f"{model_name}_{name}_fine_flux"] = flux

    if do_remove:
        delta = np.abs(
            fine_fluxes["L_final_mask"][centre]
            - fine_fluxes["A_almax_only"][centre]
        )
        target = float(compare_wave[np.argmax(delta)])
        result["remove_one_audit"] = {
            "maximum_difference_wavelength": target,
            **run_remove_one_audit(
                dll,
                sme,
                fine_wave,
                prod_ranges,
                prod_mask,
                removed_rows,
                fine_fluxes["A_almax_only"],
                fine_fluxes["L_final_mask"],
                target,
            ),
        }
    return result, arrays


def run_extra_window(model_name, linelist, label, lo, hi, workdir):
    sme, synth, setup_sec = setup_window_model(
        model_name, copy.deepcopy(linelist), lo, hi
    )
    dll = synth.get_dll()
    supported = ~np.asarray(sme.line_ion_mask, dtype=bool)
    prod_mask = np.asarray(sme.linelist["strong"], dtype=bool)
    prod_ranges = np.column_stack(
        (
            np.asarray(sme.linelist["line_range_s"], dtype=float),
            np.asarray(sme.linelist["line_range_e"], dtype=float),
        )
    )
    log = workdir / f"{model_name}-{label}-immediate.csv"
    legacy_wave, legacy_flux, legacy_sec, rows = run_prune_diagnostic(
        dll, sme, "immediate", log
    )
    removed_rows = enrich_removed_rows(rows, sme)
    removed_indices = np.asarray([row["line_index"] for row in removed_rows], dtype=int)
    final_mask = prod_mask.copy()
    final_mask[removed_indices] = False

    wfirst, wlast = float(legacy_wave[0]), float(legacy_wave[-1])
    fine_wave = np.arange(wfirst, wlast + 0.5 * FINE_STEP, FINE_STEP)
    low_almax, low_ranges = compute_low_floor_lineinfo(dll, sme, wfirst, wlast)
    valid_low = supported & np.isfinite(low_almax)
    masks = {
        "L_final_mask": final_mask,
        "A_almax_only": prod_mask,
        "P_1e-5": valid_low & (low_almax >= 1e-5),
        "F_1e-6": valid_low & (low_almax >= LOW_RANGE_FLOOR),
    }
    range_by_mode = {
        "L_final_mask": prod_ranges,
        "A_almax_only": prod_ranges,
        "P_1e-5": low_ranges,
        "F_1e-6": low_ranges,
    }
    fluxes = {}
    timings = {}
    for name in masks:
        _, flux, timing = timed_fixed_transf(
            dll, sme, fine_wave, range_by_mode[name], masks[name]
        )
        fluxes[name] = flux
        timings[name] = timing

    centre = (fine_wave >= lo) & (fine_wave <= hi)
    wave = fine_wave[centre]
    windows = {"whole_window": (lo, hi)}
    reference = fluxes["F_1e-6"][centre]
    comparisons = {
        name: metrics_on_common_grid(flux[centre], reference, wave, windows)
        for name, flux in fluxes.items()
    }
    historical = np.interp(wave, legacy_wave, legacy_flux)
    comparisons["L_historical"] = metrics_on_common_grid(
        historical, reference, wave, windows
    )
    order_result, order_fluxes = run_order_audit(
        dll, sme, workdir, fine_wave, prod_ranges, prod_mask
    )
    result = {
        "label": label,
        "window_A": [lo, hi],
        "model": model_name,
        "setup_sec": setup_sec,
        "legacy_transf_sec": legacy_sec,
        "line_counts": {
            "input_legal": int(np.sum(supported)),
            "after_almax": int(np.sum(prod_mask)),
            "almax_overlapping_transfer_window": int(
                np.sum(
                    prod_mask
                    & (prod_ranges[:, 1] > wfirst)
                    & (prod_ranges[:, 0] < wlast)
                )
            ),
            "seed_centres_tested": int(sum(row["phase"] == "seed" for row in rows)),
            "additionally_mark2": int(len(removed_rows)),
            "permissive_1e-5": int(np.sum(masks["P_1e-5"])),
            "full_1e-6": int(np.sum(masks["F_1e-6"])),
        },
        "removed_summary": {
            "by_category": count_by(removed_rows, "category"),
            "by_almax_bin": count_by(removed_rows, "almax_bin"),
            "maximum_removed_almax": float(
                max((row["almax"] for row in removed_rows), default=0.0)
            ),
            "rows": removed_rows,
        },
        "comparisons_to_F_1e-6": comparisons,
        "fine_grid_timings": timings,
        "order_audit": order_result,
    }
    arrays = {
        f"extra_{model_name}_{label}_wave": fine_wave,
        f"extra_{model_name}_{label}_legacy_wave": legacy_wave,
        f"extra_{model_name}_{label}_legacy_flux": legacy_flux,
    }
    for name, flux in fluxes.items():
        arrays[f"extra_{model_name}_{label}_{name}_flux"] = flux
    for name, flux in order_fluxes.items():
        arrays[f"extra_{model_name}_{label}_{name}_flux"] = flux
    return result, arrays


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--models", nargs="+", choices=tuple(MODELS), default=list(MODELS))
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--extra-windows", action="store_true")
    parser.add_argument("--extras-only", action="store_true")
    parser.add_argument(
        "--output", type=Path, default=ROOT / "analysis" / "rkints_second_pruning_audit.json"
    )
    parser.add_argument(
        "--npz", type=Path, default=ROOT / "analysis" / "rkints_second_pruning_audit.npz"
    )
    parser.add_argument(
        "--workdir", type=Path, default=ROOT / "analysis" / "rkints_second_pruning_logs"
    )
    args = parser.parse_args()
    args.workdir.mkdir(parents=True, exist_ok=True)
    if not BASE_BATCH_NPZ.exists():
        raise FileNotFoundError(BASE_BATCH_NPZ)
    if not LINE_LIST.exists():
        raise FileNotFoundError(LINE_LIST)

    linelist = ValdFile(str(LINE_LIST))
    with np.load(BASE_BATCH_NPZ) as source:
        base_arrays = {name: source[name] for name in source.files}
    if args.resume and args.output.exists():
        output = json.loads(args.output.read_text())
    else:
        output = {
            "setup": {
                "window_A": [LO, HI],
                "accrt": 1e-4,
                "accwi": 3e-3,
                "fine_step_A": FINE_STEP,
                "full_reference_floor": LOW_RANGE_FLOOR,
                "line_list": str(LINE_LIST),
                "diagnostic_default_unchanged": True,
            },
            "models": {},
        }
    arrays = {}
    if args.resume and args.npz.exists():
        with np.load(args.npz) as source:
            arrays.update({name: source[name] for name in source.files})

    if not args.extras_only:
        for model_name in args.models:
            print(f"running {model_name}", flush=True)
            result, model_arrays = run_model(
                model_name,
                linelist,
                args.workdir,
                base_arrays,
                do_order=model_name == "metal_poor_dwarf",
                do_remove=model_name == "metal_poor_dwarf",
            )
            output["models"][model_name] = result
            arrays.update(model_arrays)
            args.output.write_text(json.dumps(json_clean(output), indent=2, sort_keys=True) + "\n")
            np.savez_compressed(args.npz, **arrays)
            counts = result["line_counts"]
            print(
                f"{model_name}: ALMAX={counts['after_almax']} tested={counts['seed_centres_tested']} "
                f"MARK2={counts['additionally_mark2']}",
                flush=True,
            )
    if args.extra_windows:
        output.setdefault("extra_windows", {}).setdefault("metal_poor_dwarf", {})
        for label, lo, hi in (
            ("mg_b_strong", 5165.0, 5175.0),
            ("clean_atomic", 5300.0, 5310.0),
        ):
            print(f"running extra window {label}", flush=True)
            result, window_arrays = run_extra_window(
                "metal_poor_dwarf", linelist, label, lo, hi, args.workdir
            )
            output["extra_windows"]["metal_poor_dwarf"][label] = result
            arrays.update(window_arrays)
            args.output.write_text(json.dumps(json_clean(output), indent=2, sort_keys=True) + "\n")
            np.savez_compressed(args.npz, **arrays)
    print(args.output)
    print(args.npz)


if __name__ == "__main__":
    main()
