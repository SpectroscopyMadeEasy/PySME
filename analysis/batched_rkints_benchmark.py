"""Generation-batched prototype of SMElib's legacy adaptive transfer grid.

This diagnostic intentionally does not change the production RKINTS path.
It constructs the legacy line-centre seed grid in Python and evaluates each
refinement generation with one fixed-grid ``Transf`` call.  The prototype
uses the same midpoint error criterion and 0.3 km/s minimum spacing as
RKINTS, while making the unavoidable line-state difference explicit: legacy
RKINTS can deactivate a line immediately after evaluating its centre,
whereas a generation batch keeps the precomputed ALMAX mask fixed.
"""

from __future__ import annotations

import argparse
import copy
import json
import time
from pathlib import Path

import numpy as np

from almax_seeded_transfer_grid_benchmark import (
    C_KMS,
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
    post_convolve,
)
from pysme.linelist.vald import ValdFile
from pysme.synthesize import Synthesizer


DVEL_MIN_KMS = 0.3
COMMON_STEP = 0.000625
NWSIZE = 40000


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


def raw_transf(dll, sme, wave, *, keep_lineop, adaptive_mode=None):
    if adaptive_mode is not None:
        dll.SetAdaptiveTransferGridMode(adaptive_mode)
    started = time.perf_counter()
    _, got_wave, sint, cint = dll.Transf(
        sme.mu,
        wave=wave,
        nwmax=NWSIZE if wave is None else len(wave),
        accrt=sme.accrt,
        accwi=sme.accwi,
        keep_lineop=keep_lineop,
        long_continuum=True,
    )
    return (
        np.asarray(got_wave, dtype=float),
        np.asarray(sint, dtype=float),
        np.asarray(cint, dtype=float),
        time.perf_counter() - started,
    )


def normalized_intrinsic(sme, sint, cint):
    line = integrate_mu_areas(sme.mu, sint)
    continuum = integrate_mu_areas(sme.mu, cint)
    return np.asarray(line / continuum, dtype=float)


def restore_precomputed_line_state(dll, sme):
    """Restore the strict ALMAX mask/ranges after legacy mutates MARK."""
    supported = ~np.asarray(sme.line_ion_mask, dtype=bool)
    ranges_s = np.asarray(sme.linelist["line_range_s"], dtype=float)[supported]
    ranges_e = np.asarray(sme.linelist["line_range_e"], dtype=float)[supported]
    strong = np.asarray(sme.linelist["strong"], dtype=np.uint8)[supported]
    dll.SetLineInfoMode(2)
    dll.InputLinePrecomputedInfo(ranges_s, ranges_e, strong)


def legacy_initial_grid(wfirst, wlast, line_wave, active):
    """Generate RKINTS' initial endpoints/midpoints/line-centre sequence."""
    nodes = [float(wfirst)]
    accepted_centres = 0
    for centre, enabled in zip(line_wave, active):
        centre = float(centre)
        if not enabled:
            continue
        min_step = centre * DVEL_MIN_KMS / C_KMS
        if centre > wfirst and centre < wlast and centre - nodes[-1] > min_step:
            nodes.append(0.5 * (centre + nodes[-1]))
            nodes.append(centre)
            accepted_centres += 1

    min_step = float(wlast) * DVEL_MIN_KMS / C_KMS
    if float(wlast) - nodes[-1] > min_step:
        nodes.append(float(wlast))
    else:
        # This is the slightly surprising overwrite used by RKINTS.
        nodes[-1] = float(wlast)
    return np.asarray(nodes, dtype=float), accepted_centres


def _cache_columns(cache, wave, sint, cint):
    for column, wavelength in enumerate(wave):
        cache[float(wavelength)] = (
            np.asarray(sint[:, column], dtype=float).copy(),
            np.asarray(cint[:, column], dtype=float).copy(),
        )


def batched_refinement(dll, sme, initial_wave):
    """Evaluate all active interval midpoints once per refinement level."""
    cache = {}
    generations = []
    bookkeeping_sec = 0.0

    wave, sint, cint, elapsed = raw_transf(
        dll, sme, np.asarray(initial_wave, dtype=float), keep_lineop=False
    )
    _cache_columns(cache, wave, sint, cint)
    generations.append(
        {
            "generation": 0,
            "kind": "initial_grid",
            "requested": int(initial_wave.size),
            "evaluated": int(wave.size),
            "cache_hits": 0,
            "transf_sec": elapsed,
            "active_intervals_after": int(max(0, wave.size - 1)),
        }
    )

    active_intervals = [(float(left), float(right)) for left, right in zip(wave[:-1], wave[1:])]
    imu_ref = int(np.argmax(np.asarray(sme.mu, dtype=float)))
    accepted_intervals = 0
    unresolved_intervals = 0
    generation_number = 1

    while active_intervals:
        bookkeeping_started = time.perf_counter()
        midpoints = [0.5 * (left + right) for left, right in active_intervals]
        unique_midpoints = np.asarray(sorted(set(midpoints)), dtype=float)
        missing = np.asarray(
            [value for value in unique_midpoints if float(value) not in cache],
            dtype=float,
        )
        bookkeeping_sec += time.perf_counter() - bookkeeping_started

        if len(cache) + missing.size >= NWSIZE:
            unresolved_intervals += len(active_intervals)
            break

        transfer_sec = 0.0
        evaluated = 0
        if missing.size:
            got_wave, got_sint, got_cint, transfer_sec = raw_transf(
                dll, sme, missing, keep_lineop=True
            )
            _cache_columns(cache, got_wave, got_sint, got_cint)
            evaluated = int(got_wave.size)

        bookkeeping_started = time.perf_counter()
        next_intervals = []
        generation_worst = 0.0
        for left, right in active_intervals:
            midpoint = 0.5 * (left + right)
            line_left, _ = cache[left]
            line_mid, cont_mid = cache[midpoint]
            line_right, _ = cache[right]
            fnorm = float(cont_mid[imu_ref])
            numerator = abs(
                float(line_mid[imu_ref])
                - 0.5 * (float(line_left[imu_ref]) + float(line_right[imu_ref]))
            ) + 0.005 * abs(float(line_left[imu_ref]) - float(line_right[imu_ref]))
            error = numerator / max(abs(fnorm), np.finfo(float).tiny)
            generation_worst = max(generation_worst, error)
            # RKINTS defines DWL_MIN from the newly inserted midpoint.
            minimum_half_width = midpoint * DVEL_MIN_KMS / C_KMS
            if error < sme.accwi or midpoint - left <= minimum_half_width:
                accepted_intervals += 1
            else:
                next_intervals.append((left, midpoint))
                next_intervals.append((midpoint, right))

        bookkeeping_sec += time.perf_counter() - bookkeeping_started
        generations.append(
            {
                "generation": generation_number,
                "kind": "midpoint_probes",
                "requested": int(len(midpoints)),
                "evaluated": evaluated,
                "cache_hits": int(len(unique_midpoints) - evaluated),
                "transf_sec": transfer_sec,
                "worst_error": generation_worst,
                "accepted_intervals": int(len(active_intervals) - len(next_intervals) // 2),
                "active_intervals_after": int(len(next_intervals)),
            }
        )
        active_intervals = next_intervals
        generation_number += 1

    final_wave = np.asarray(sorted(cache), dtype=float)
    final_sint = np.column_stack([cache[value][0] for value in final_wave])
    final_cint = np.column_stack([cache[value][1] for value in final_wave])
    return final_wave, final_sint, final_cint, {
        "generation_count_including_initial": int(len(generations)),
        "midpoint_generations": int(max(0, len(generations) - 1)),
        "native_transf_calls": int(len(generations)),
        "native_wavelength_evaluations": int(sum(item["evaluated"] for item in generations)),
        "transf_sec": float(sum(item["transf_sec"] for item in generations)),
        "bookkeeping_sec": float(bookkeeping_sec),
        "total_sec": float(sum(item["transf_sec"] for item in generations) + bookkeeping_sec),
        "accepted_intervals": int(accepted_intervals),
        "unresolved_intervals": int(unresolved_intervals),
        "generations": generations,
    }


def nearest_grid_distance(source, target):
    """Distance from every source node to its nearest target node."""
    source = np.asarray(source, dtype=float)
    target = np.asarray(target, dtype=float)
    position = np.searchsorted(target, source)
    left = target[np.clip(position - 1, 0, target.size - 1)]
    right = target[np.clip(position, 0, target.size - 1)]
    return np.minimum(np.abs(source - left), np.abs(source - right))


def match_common_nodes(left, right, atol=5e-11):
    left = np.asarray(left, dtype=float)
    right = np.asarray(right, dtype=float)
    left_index = []
    right_index = []
    i = j = 0
    while i < left.size and j < right.size:
        delta = left[i] - right[j]
        if abs(delta) <= atol:
            left_index.append(i)
            right_index.append(j)
            i += 1
            j += 1
        elif delta < 0:
            i += 1
        else:
            j += 1
    return np.asarray(left_index, dtype=int), np.asarray(right_index, dtype=int)


def grid_comparison(legacy_wave, batch_wave):
    li, bi = match_common_nodes(legacy_wave, batch_wave)
    legacy_distance = nearest_grid_distance(legacy_wave, batch_wave)
    batch_distance = nearest_grid_distance(batch_wave, legacy_wave)
    tolerance = 5e-11
    return {
        "legacy_points": int(legacy_wave.size),
        "batch_points": int(batch_wave.size),
        "common_points": int(li.size),
        "legacy_only_points": int(np.sum(legacy_distance > tolerance)),
        "batch_only_points": int(np.sum(batch_distance > tolerance)),
        "max_legacy_to_batch_distance_A": float(np.max(legacy_distance)),
        "max_batch_to_legacy_distance_A": float(np.max(batch_distance)),
        "hausdorff_distance_A": float(max(np.max(legacy_distance), np.max(batch_distance))),
        "common_legacy_indices": li,
        "common_batch_indices": bi,
    }


def spectrum_metrics(wave, flux, reference_wave, reference_flux, windows):
    candidate = np.interp(reference_wave, wave, flux)
    result = {
        "intrinsic": error_stats(candidate, reference_flux),
        "equivalent_width_A": equivalent_widths(reference_wave, candidate, windows),
        "post_convolution": {},
    }
    step = float(reference_wave[1] - reference_wave[0])
    for resolving_power in (20000.0, 60000.0):
        ref_conv = post_convolve(reference_flux, step, resolving_power)
        got_conv = post_convolve(candidate, step, resolving_power)
        result["post_convolution"][str(int(resolving_power))] = error_stats(got_conv, ref_conv)
    return result, candidate


def paired_spectrum_metrics(wave, candidate, reference, windows):
    """Compare two spectra already sampled on the same regular grid."""
    result = {
        "intrinsic": error_stats(candidate, reference),
        "equivalent_width_delta_A": {},
        "post_convolution": {},
    }
    candidate_ew = equivalent_widths(wave, candidate, windows)
    reference_ew = equivalent_widths(wave, reference, windows)
    result["equivalent_width_delta_A"] = {
        name: candidate_ew[name] - reference_ew[name] for name in windows
    }
    step = float(wave[1] - wave[0])
    for resolving_power in (20000.0, 60000.0):
        result["post_convolution"][str(int(resolving_power))] = error_stats(
            post_convolve(candidate, step, resolving_power),
            post_convolve(reference, step, resolving_power),
        )
    return result


def run_model(model_name, linelist):
    setup_wave = np.arange(LO, HI + 0.025, 0.05)
    sme = make_sme(model_name, linelist, setup_wave)
    synth = Synthesizer()
    setup_started = time.perf_counter()
    synth.synthesize_spectrum(
        sme,
        updateStructure=False,
        reuse_wavelength_grid=True,
        radial_velocity_mode="fast",
        linelist_mode="all",
        smelib_lineinfo_mode=2,
    )
    setup_sec = time.perf_counter() - setup_started
    dll = synth.get_dll()
    restore_precomputed_line_state(dll, sme)

    # Current production default: strict ALMAX information plus legacy RKINTS.
    legacy_wave, legacy_sint, legacy_cint, legacy_sec = raw_transf(
        dll, sme, None, keep_lineop=False, adaptive_mode="legacy"
    )
    legacy_flux = normalized_intrinsic(sme, legacy_sint, legacy_cint)

    line_wave = np.asarray(sme.linelist["wlcent"], dtype=float)
    supported = ~np.asarray(sme.line_ion_mask, dtype=bool)
    strong = np.asarray(sme.linelist["strong"], dtype=bool)
    active = supported & strong
    initial_wave, accepted_centres = legacy_initial_grid(
        float(legacy_wave[0]), float(legacy_wave[-1]), line_wave, active
    )

    restore_precomputed_line_state(dll, sme)
    batch_wave, batch_sint, batch_cint, batch_timing = batched_refinement(
        dll, sme, initial_wave
    )
    batch_flux = normalized_intrinsic(sme, batch_sint, batch_cint)

    comparison = grid_comparison(legacy_wave, batch_wave)
    legacy_common = comparison.pop("common_legacy_indices")
    batch_common = comparison.pop("common_batch_indices")
    common_flux_error = error_stats(
        batch_flux[batch_common], legacy_flux[legacy_common]
    )

    common_wave = np.arange(LO, HI + 0.5 * COMMON_STEP, COMMON_STEP)
    # The fixed-grid reference shares the batch's strict, immutable ALMAX mask.
    reference_wave, reference_sint, reference_cint, reference_sec = raw_transf(
        dll, sme, common_wave, keep_lineop=True
    )
    reference_flux = normalized_intrinsic(sme, reference_sint, reference_cint)

    almax = np.asarray(sme.linelist["almax_ratio"], dtype=float)
    molecule = is_molecular(np.asarray(sme.linelist["species"]))
    windows = feature_windows(line_wave, almax, molecule)
    reference_ew = equivalent_widths(reference_wave, reference_flux, windows)
    legacy_metrics, legacy_resampled = spectrum_metrics(
        legacy_wave, legacy_flux, reference_wave, reference_flux, windows
    )
    batch_metrics, batch_resampled = spectrum_metrics(
        batch_wave, batch_flux, reference_wave, reference_flux, windows
    )
    legacy_metrics["equivalent_width_delta_A"] = {
        name: legacy_metrics["equivalent_width_A"][name] - reference_ew[name]
        for name in windows
    }
    batch_metrics["equivalent_width_delta_A"] = {
        name: batch_metrics["equivalent_width_A"][name] - reference_ew[name]
        for name in windows
    }
    batch_vs_legacy = paired_spectrum_metrics(
        reference_wave, batch_resampled, legacy_resampled, windows
    )

    result = {
        "model": model_name,
        "stellar_parameters": dict(zip(("teff", "logg", "monh", "vmic"), MODELS[model_name])),
        "setup_sec": setup_sec,
        "line_state": {
            "input_lines": int(len(sme.linelist)),
            "supported_lines": int(np.sum(supported)),
            "precomputed_strong_lines": int(np.sum(active)),
            "accepted_seed_centres": int(accepted_centres),
            "sequential_dependency": (
                "legacy RKINTS can set MARK=2 after each line centre; the batch keeps "
                "the strict precomputed ALMAX mask fixed for every generation"
            ),
        },
        "legacy": {
            "transfer_points": int(legacy_wave.size),
            "transf_sec": legacy_sec,
            "native_transf_calls": 1,
            "conceptual_single_wavelength_rt_evaluations": int(legacy_wave.size),
            "metrics_vs_fixed_reference": legacy_metrics,
        },
        "batched": {
            "initial_points": int(initial_wave.size),
            "transfer_points": int(batch_wave.size),
            **batch_timing,
            "speedup_vs_legacy": float(legacy_sec / batch_timing["total_sec"]),
            "metrics_vs_fixed_reference": batch_metrics,
        },
        "fixed_reference": {
            "spacing_A": COMMON_STEP,
            "transfer_points": int(reference_wave.size),
            "transf_sec": reference_sec,
            "equivalent_width_A": reference_ew,
        },
        "grid_comparison": comparison,
        "common_node_batch_minus_legacy": common_flux_error,
        "batch_minus_legacy": batch_vs_legacy,
        "feature_windows": {name: list(bounds) for name, bounds in windows.items()},
    }
    arrays = {
        f"{model_name}_legacy_wave": legacy_wave,
        f"{model_name}_legacy_flux": legacy_flux,
        f"{model_name}_batch_wave": batch_wave,
        f"{model_name}_batch_flux": batch_flux,
        f"{model_name}_initial_wave": initial_wave,
        f"{model_name}_reference_wave": reference_wave,
        f"{model_name}_reference_flux": reference_flux,
        f"{model_name}_legacy_resampled": legacy_resampled,
        f"{model_name}_batch_resampled": batch_resampled,
    }
    return result, arrays


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--models", nargs="+", choices=tuple(MODELS), default=list(MODELS))
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("analysis/batched_rkints_benchmark.json"),
    )
    parser.add_argument(
        "--npz",
        type=Path,
        default=Path("analysis/batched_rkints_benchmark.npz"),
    )
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()
    if not LINE_LIST.exists():
        raise FileNotFoundError(LINE_LIST)

    linelist = ValdFile(str(LINE_LIST))
    if args.resume and args.output.exists():
        output = json.loads(args.output.read_text())
    else:
        output = {
            "setup": {
                "window_A": [LO, HI],
                "common_reference_step_A": COMMON_STEP,
                "accrt": 1e-4,
                "accwi": 3e-3,
                "minimum_velocity_step_kms": DVEL_MIN_KMS,
                "line_list": str(LINE_LIST),
                "prototype": "Python generation scheduler + existing optimized fixed-grid Transf",
            },
            "models": {},
        }
    arrays = {}
    if args.resume and args.npz.exists():
        with np.load(args.npz) as existing:
            arrays.update({name: existing[name] for name in existing.files})

    for model_name in args.models:
        print(f"running {model_name}", flush=True)
        result, model_arrays = run_model(model_name, copy.deepcopy(linelist))
        output["models"][model_name] = result
        arrays.update(model_arrays)
        args.output.write_text(json.dumps(json_clean(output), indent=2, sort_keys=True) + "\n")
        np.savez_compressed(args.npz, **arrays)
        print(
            f"{model_name}: legacy={result['legacy']['transf_sec']:.3f}s "
            f"batch={result['batched']['total_sec']:.3f}s "
            f"speedup={result['batched']['speedup_vs_legacy']:.2f}x",
            flush=True,
        )

    print(args.output)
    print(args.npz)


if __name__ == "__main__":
    main()
