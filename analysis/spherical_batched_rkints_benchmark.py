"""Audit and benchmark native generation-batched spherical RKINTS.

The benchmark uses real spherical MARCS atmospheres and the same immutable
ALMAX mask/physical ranges for the controlled batched and fixed-grid paths.
Historical legacy RKINTS_sph is run separately because it may shorten Wlim
during its line-centre seed scan.
"""

from __future__ import annotations

import argparse
import ctypes
import json
import os
import re
import tempfile
import time
from pathlib import Path

import numpy as np
from scipy.ndimage import gaussian_filter1d

from line_centric_wide_benchmark import extract_vald_window
from pysme.atmosphere.marcsfile import MarcsAtmosphere
from pysme.linelist.vald import ValdFile
from pysme.sme import SME_Structure
from pysme.sme_synth import SME_DLL
from pysme.synthesize import Synthesizer


ROOT = Path(__file__).resolve().parents[1]
SOURCE_LINES = Path("/Users/mingjie/Documents/Research/data/merged_3700-9500_hfs.lin")
DENSE_LINES = Path("/tmp/pysme_line_centric_windows/window_0010A_pad_150A.lin")
CACHE = Path("/tmp/pysme_spherical_batched")
MODEL_ROOT = Path("/Users/mingjie/code/TSFitPy/input_files/model_atmospheres/1D")
C_KMS = 299792.458
DVEL_MIN_KMS = 0.3
ACCRT = 1e-4
RANGE_FLOOR = 1e-6
ACCWI = 3e-3
MU = np.polynomial.legendre.leggauss(7)[0]
MU = 0.5 * (MU + 1.0)

MODELS = {
    "moderate_giant": MODEL_ROOT / "s4500_g+2.0_m1.0_t02_st_z+0.00_a+0.00_c+0.00_n+0.00_o+0.00_r+0.00_s+0.00.mod",
    "cool_line_rich_giant": MODEL_ROOT / "s3500_g+1.0_m1.0_t02_st_z+0.00_a+0.00_c+0.00_n+0.00_o+0.00_r+0.00_s+0.00.mod",
    "metal_poor_giant": MODEL_ROOT / "s4500_g+2.0_m1.0_t02_st_z-2.00_a+0.40_c+0.00_n+0.00_o+0.40_r+0.00_s+0.00.mod",
}

CASES = {
    "dense_moderate": ("moderate_giant", 5195.0, 5205.0, "dense"),
    "dense_cool": ("cool_line_rich_giant", 5195.0, 5205.0, "dense"),
    "dense_metal_poor": ("metal_poor_giant", 5195.0, 5205.0, "dense"),
    "mgb_moderate": ("moderate_giant", 5165.0, 5190.0, "dense"),
    "halpha_moderate": ("moderate_giant", 6550.0, 6575.0, "halpha"),
}


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


def ensure_halpha_lines() -> Path:
    CACHE.mkdir(parents=True, exist_ok=True)
    destination = CACHE / "halpha_6250_6875.lin"
    if not destination.exists():
        extract_vald_window(SOURCE_LINES, destination, 6250.0, 6875.0)
    return destination


def parse_timing(text: str) -> dict:
    starts = [match.start() for match in re.finditer(r"SMElib timing \[Transf#", text)]
    block = text[starts[-1] :] if starts else text
    result = {"raw": block.strip()}
    scalar = {
        "total_sec": r"^  total\s+:\s+([0-9.eE+-]+) s",
        "setup_sec": r"^  setup\s+:\s+([0-9.eE+-]+) s",
        "lineopac_sec": r"^  LINEOPAC\s+:\s+([0-9.eE+-]+) s, calls=(\d+)",
        "range_scan_sec": r"^  range_scan\s+:\s+([0-9.eE+-]+) s, scans=(\d+), lines=(\d+)",
        "rkints_sph_sec": r"^  RKINTS_sph\s+:\s+([0-9.eE+-]+) s, calls=(\d+)",
        "opmtrx_sec": r"^  OPMTRX\s+:\s+([0-9.eE+-]+) s, calls=(\d+)",
        "tbintg_sph_sec": r"^  TBINTG_sph\s+:\s+([0-9.eE+-]+) s, calls=(\d+)",
    }
    for name, pattern in scalar.items():
        match = re.search(pattern, block, flags=re.MULTILINE)
        if match:
            result[name] = float(match.group(1))
            if match.lastindex and match.lastindex >= 2:
                result[name.replace("_sec", "_calls")] = int(match.group(2))
            if name == "range_scan_sec" and match.lastindex == 3:
                result["range_scan_lines"] = int(match.group(3))
    match = re.search(
        r"^  intervals\s+: calls=(\d+), candidates=(\d+), full_scan=(\d+), reduction=([0-9.eE+-]+)%",
        block,
        flags=re.MULTILINE,
    )
    if match:
        result.update(
            interval_calls=int(match.group(1)),
            interval_candidates=int(match.group(2)),
            interval_full_scan=int(match.group(3)),
            interval_reduction_fraction=float(match.group(4)) / 100.0,
        )
    match = re.search(
        r"^  adaptive\s+: generations=(\d+), probes=(\d+), refined=(\d+), accepted=(\d+)",
        block,
        flags=re.MULTILINE,
    )
    if match:
        result.update(
            adaptive_generations=int(match.group(1)),
            adaptive_probes=int(match.group(2)),
            adaptive_intervals_refined=int(match.group(3)),
            adaptive_intervals_accepted=int(match.group(4)),
        )
    return result


def capture_transf(dll, *, mode, wave=None):
    dll.SetAdaptiveTransferGridMode(mode)
    with tempfile.TemporaryFile(mode="w+b") as stream:
        saved_stderr = os.dup(2)
        try:
            os.dup2(stream.fileno(), 2)
            started = time.perf_counter()
            nw, got_wave, sint, cint = dll.Transf(
                MU,
                wave=wave,
                nwmax=400000 if wave is None else len(wave),
                accrt=ACCRT,
                accwi=ACCWI,
                keep_lineop=False,
                long_continuum=True,
            )
            wall_sec = time.perf_counter() - started
            ctypes.CDLL(None).fflush(None)
        finally:
            os.dup2(saved_stderr, 2)
            os.close(saved_stderr)
        stream.seek(0)
        timing_text = stream.read().decode("utf-8", errors="replace")
    return (
        int(nw),
        np.asarray(got_wave, dtype=float),
        np.asarray(sint, dtype=float),
        np.asarray(cint, dtype=float),
        wall_sec,
        parse_timing(timing_text),
    )


def prepare_state(dll, atmosphere, linelist, lo, hi):
    dll.SetLibraryPath()
    dll.SetContinuumOpacityGrid("adaptive", rtol=1e-3)
    ion_mask = np.asarray(dll.InputLineList(linelist), dtype=bool)
    dll.InputModel(atmosphere.teff, atmosphere.logg, atmosphere.vturb, atmosphere)
    dll.InputAbund(atmosphere.abund)
    dll.Ionization(0)
    dll.SetVWscale(1.0)
    dll.SetH2broad()
    dll.InputWaveRange(lo, hi)
    dll.Opacity()
    almax, ranges = dll.ALMAXRange(accrt=RANGE_FLOOR)
    almax = np.asarray(almax, dtype=float)
    ranges = np.asarray(ranges, dtype=float)
    strong = np.asarray(almax >= ACCRT, dtype=np.uint8)
    dll.InputLinePrecomputedInfo(ranges[:, 0], ranges[:, 1], strong)
    dll.SetLineInfoMode(2)
    return {
        "ion_mask": ion_mask,
        "almax": almax,
        "ranges": ranges,
        "strong": strong,
        "supported_wave": np.asarray(linelist["wlcent"], dtype=float)[~ion_mask],
    }


def run_mode(dll, atmosphere, linelist, lo, hi, mode, wave=None):
    state = prepare_state(dll, atmosphere, linelist, lo, hi)
    before = state["ranges"].copy()
    nw, got_wave, sint, cint, wall, timing = capture_transf(
        dll, mode=mode, wave=wave
    )
    after = np.asarray(dll.GetLineRange(), dtype=float)
    return {
        "nw": nw,
        "wave": got_wave,
        "sint": sint,
        "cint": cint,
        "wall_sec": wall,
        "timing": timing,
        "ranges_before": before,
        "ranges_after": after,
        "state": state,
    }


def initial_grid(lo, hi, line_wave, active):
    nodes = [float(lo)]
    for centre, enabled in zip(line_wave, active):
        centre = float(centre)
        if (
            enabled
            and lo < centre < hi
            and centre - nodes[-1] > centre * DVEL_MIN_KMS / C_KMS
        ):
            nodes.append(0.5 * (centre + nodes[-1]))
            nodes.append(centre)
    if hi - nodes[-1] > hi * DVEL_MIN_KMS / C_KMS:
        nodes.append(float(hi))
    else:
        nodes[-1] = float(hi)
    return np.asarray(nodes)


def replay_sequential_grid(batch, initial):
    """Replay depth-first sequential decisions from immutable node values."""
    wave = batch["wave"]
    sint = batch["sint"]
    cint = batch["cint"]
    imu_ref = int(np.argmax(MU))
    used = set(float(value) for value in initial)

    def node_index(value):
        index = int(np.searchsorted(wave, value))
        candidates = [i for i in (index - 1, index) if 0 <= i < wave.size]
        best = min(candidates, key=lambda i: abs(wave[i] - value))
        if abs(wave[best] - value) > 5e-11:
            raise AssertionError(f"missing controlled midpoint {value}")
        return best

    def refine(left, right):
        midpoint = 0.5 * (left + right)
        li, mi, ri = node_index(left), node_index(midpoint), node_index(right)
        used.add(float(wave[mi]))
        continuum = max(abs(float(cint[imu_ref, mi])), np.finfo(float).tiny)
        error = (
            abs(float(sint[imu_ref, mi]) - 0.5 * (sint[imu_ref, li] + sint[imu_ref, ri]))
            + 0.005 * abs(float(sint[imu_ref, li] - sint[imu_ref, ri]))
        ) / continuum
        if error >= ACCWI and midpoint - left > left * DVEL_MIN_KMS / C_KMS:
            refine(left, midpoint)
            refine(midpoint, right)

    for left, right in zip(initial[:-1], initial[1:]):
        refine(float(left), float(right))
    return np.asarray(sorted(used), dtype=float)


def nearest_distance(source, target):
    position = np.searchsorted(target, source)
    left = target[np.clip(position - 1, 0, target.size - 1)]
    right = target[np.clip(position, 0, target.size - 1)]
    return np.minimum(np.abs(source - left), np.abs(source - right))


def grid_metrics(legacy, batch):
    tolerance = 5e-11
    ld = nearest_distance(legacy, batch)
    bd = nearest_distance(batch, legacy)
    return {
        "legacy_points": int(legacy.size),
        "batch_points": int(batch.size),
        "common_points": int(np.sum(ld <= tolerance)),
        "legacy_only": int(np.sum(ld > tolerance)),
        "batch_only": int(np.sum(bd > tolerance)),
        "hausdorff_A": float(max(np.max(ld), np.max(bd))),
    }


def common_node_errors(left, right):
    li, ri = [], []
    i = j = 0
    while i < left["wave"].size and j < right["wave"].size:
        delta = left["wave"][i] - right["wave"][j]
        if abs(delta) <= 5e-11:
            li.append(i)
            ri.append(j)
            i += 1
            j += 1
        elif delta < 0:
            i += 1
        else:
            j += 1
    li, ri = np.asarray(li), np.asarray(ri)
    ds = right["sint"][:, ri] - left["sint"][:, li]
    dc = right["cint"][:, ri] - left["cint"][:, li]
    return {
        "points": int(li.size),
        "sint_max_abs": float(np.max(np.abs(ds))),
        "sint_rms": float(np.sqrt(np.mean(ds * ds))),
        "cint_max_abs": float(np.max(np.abs(dc))),
        "cint_rms": float(np.sqrt(np.mean(dc * dc))),
        "per_ray_sint_max_abs": np.max(np.abs(ds), axis=1),
        "per_ray_cint_max_abs": np.max(np.abs(dc), axis=1),
    }


def integrate_mu(intensity):
    radius = np.sqrt(1.0 - MU * MU)
    order = np.argsort(radius)
    radius = radius[order]
    values = intensity[order]
    boundaries = np.sqrt(0.5 * (radius[:-1] ** 2 + radius[1:] ** 2))
    boundaries = np.concatenate(([0.0], boundaries, [1.0]))
    weights = boundaries[1:] ** 2 - boundaries[:-1] ** 2
    return np.pi * np.sum(weights[:, None] * values, axis=0)


def regular_flux(result, lo, hi, dv=0.3):
    log_step = dv / C_KMS
    log_wave = np.arange(np.log(lo), np.log(hi) + 0.5 * log_step, log_step)
    wave = np.exp(log_wave)
    sint = np.vstack([np.interp(wave, result["wave"], row) for row in result["sint"]])
    cint = np.vstack([np.interp(wave, result["wave"], row) for row in result["cint"]])
    return wave, integrate_mu(sint) / integrate_mu(cint)


def flux_metrics(legacy, batch, lo, hi):
    wave, legacy_flux = regular_flux(legacy, lo, hi)
    _, batch_flux = regular_flux(batch, lo, hi)
    delta = batch_flux - legacy_flux
    result = {
        "regular_velocity_step_kms": 0.3,
        "intrinsic_max_abs": float(np.max(np.abs(delta))),
        "intrinsic_rms": float(np.sqrt(np.mean(delta * delta))),
        "legacy_ew_A": float(np.trapezoid(1.0 - legacy_flux, wave)),
        "batch_ew_A": float(np.trapezoid(1.0 - batch_flux, wave)),
        "ew_delta_A": float(np.trapezoid(legacy_flux - batch_flux, wave)),
        "convolved": {},
    }
    for resolution in (20000, 60000):
        sigma_pixels = (C_KMS / resolution / 2.354820045) / 0.3
        lconv = gaussian_filter1d(legacy_flux, sigma_pixels, mode="nearest")
        bconv = gaussian_filter1d(batch_flux, sigma_pixels, mode="nearest")
        dconv = bconv - lconv
        result["convolved"][str(resolution)] = {
            "max_abs": float(np.max(np.abs(dconv))),
            "rms": float(np.sqrt(np.mean(dconv * dconv))),
        }
    return result, wave, legacy_flux, batch_flux


def run_case(name, linelists):
    model_name, lo, hi, line_key = CASES[name]
    atmosphere = MarcsAtmosphere(str(MODELS[model_name]))
    linelist = linelists[line_key]
    dll = SME_DLL()

    legacy = run_mode(dll, atmosphere, linelist, lo, hi, "legacy")
    batch = run_mode(dll, atmosphere, linelist, lo, hi, "batched")
    fixed = run_mode(dll, atmosphere, linelist, lo, hi, "batched", wave=batch["wave"])

    initial = initial_grid(
        lo, hi, batch["state"]["supported_wave"], batch["state"]["strong"]
    )
    sequential_wave = replay_sequential_grid(batch, initial)
    controlled_sint = fixed["sint"] - batch["sint"]
    controlled_cint = fixed["cint"] - batch["cint"]
    flux, regular_wave, legacy_flux, batch_flux = flux_metrics(legacy, batch, lo, hi)

    timing_total = batch["timing"].get("total_sec", batch["wall_sec"])
    legacy_total = legacy["timing"].get("total_sec", legacy["wall_sec"])
    result = {
        "case": name,
        "model": model_name,
        "atmosphere": {
            "path": str(MODELS[model_name]),
            "teff": atmosphere.teff,
            "logg": atmosphere.logg,
            "monh": atmosphere.monh,
            "geometry": atmosphere.geom,
            "layers": len(atmosphere.rhox),
            "radius_cm": atmosphere.radius,
        },
        "window_A": [lo, hi],
        "input_lines": len(linelist),
        "supported_lines": int(batch["state"]["almax"].size),
        "active_lines": int(np.sum(batch["state"]["strong"])),
        "mu": MU,
        "legacy": {
            "points": legacy["nw"],
            "wall_sec": legacy["wall_sec"],
            "timing": legacy["timing"],
            "active_range_values_changed": int(
                np.sum(
                    legacy["ranges_before"][legacy["state"]["strong"].astype(bool)]
                    != legacy["ranges_after"][legacy["state"]["strong"].astype(bool)]
                )
            ),
        },
        "batch": {
            "points": batch["nw"],
            "wall_sec": batch["wall_sec"],
            "timing": batch["timing"],
            "active_range_values_changed": int(
                np.sum(
                    batch["ranges_before"][batch["state"]["strong"].astype(bool)]
                    != batch["ranges_after"][batch["state"]["strong"].astype(bool)]
                )
            ),
            "speedup": float(legacy_total / timing_total),
        },
        "grid_historical": grid_metrics(legacy["wave"], batch["wave"]),
        "historical_common_nodes": common_node_errors(legacy, batch),
        "controlled_semantics": {
            "sequential_replay_points": int(sequential_wave.size),
            "grid_bitwise_equal": bool(np.array_equal(sequential_wave, batch["wave"])),
            "fixed_common_wave_bitwise_equal": bool(
                np.array_equal(fixed["wave"], batch["wave"])
            ),
            "sint_bitwise_equal": bool(np.array_equal(fixed["sint"], batch["sint"])),
            "cint_bitwise_equal": bool(np.array_equal(fixed["cint"], batch["cint"])),
            "sint_max_abs": float(np.max(np.abs(controlled_sint))),
            "cint_max_abs": float(np.max(np.abs(controlled_cint))),
            "per_ray_sint_max_abs": np.max(np.abs(controlled_sint), axis=1),
            "per_ray_cint_max_abs": np.max(np.abs(controlled_cint), axis=1),
        },
        "flux_batch_minus_historical_legacy": flux,
    }
    arrays = {
        f"{name}_legacy_wave": legacy["wave"],
        f"{name}_legacy_sint": legacy["sint"],
        f"{name}_legacy_cint": legacy["cint"],
        f"{name}_batch_wave": batch["wave"],
        f"{name}_batch_sint": batch["sint"],
        f"{name}_batch_cint": batch["cint"],
        f"{name}_sequential_replay_wave": sequential_wave,
        f"{name}_regular_wave": regular_wave,
        f"{name}_legacy_flux": legacy_flux,
        f"{name}_batch_flux": batch_flux,
    }
    dll.SetLineInfoMode(0)
    dll.SetAdaptiveTransferGridMode("batched")
    return result, arrays


def run_spherical_halpha_nlte():
    """Use the cached, tested H departure grid with a real spherical MARCS model."""
    atmosphere = MarcsAtmosphere(str(MODELS["moderate_giant"]))
    sme = SME_Structure()
    sme.teff = atmosphere.teff
    sme.logg = atmosphere.logg
    sme.monh = atmosphere.monh
    sme.vmic = atmosphere.vturb
    sme.vmac = 0.0
    sme.vsini = 0.0
    sme.abund = atmosphere.abund
    sme.atmo = atmosphere
    sme.linelist = ValdFile(str(ROOT / "test/halpha_window_cdr_union.lin"))
    sme.wran = [[6561.0, 6564.2]]
    sme.vrad_flag = "none"
    sme.cscale_flag = "none"
    sme.accrt = ACCRT
    sme.accwi = ACCWI
    sme.transfer_grid_method = "batched"
    sme.nlte.set_nlte("H", "nlte_H_pysme.grd")

    synthesizer = Synthesizer()
    synthesizer.synthesize_spectrum(
        sme,
        updateStructure=False,
        reuse_wavelength_grid=True,
        radial_velocity_mode="fast",
        linelist_mode="all",
        smelib_lineinfo_mode=2,
    )
    dll = synthesizer.get_dll()
    supported = ~np.asarray(sme.line_ion_mask, dtype=bool)
    range_s = np.asarray(sme.linelist["line_range_s"], dtype=float)[supported]
    range_e = np.asarray(sme.linelist["line_range_e"], dtype=float)[supported]
    strong = np.asarray(sme.linelist["strong"], dtype=np.uint8)[supported]

    def transfer(mode, wave=None):
        dll.InputLinePrecomputedInfo(range_s, range_e, strong)
        dll.SetLineInfoMode(2)
        return capture_transf(dll, mode=mode, wave=wave)

    legacy_values = transfer("legacy")
    batch_values = transfer("batched")
    fixed_values = transfer("batched", wave=batch_values[1])
    legacy = {
        "nw": legacy_values[0], "wave": legacy_values[1],
        "sint": legacy_values[2], "cint": legacy_values[3],
        "wall_sec": legacy_values[4], "timing": legacy_values[5],
    }
    batch = {
        "nw": batch_values[0], "wave": batch_values[1],
        "sint": batch_values[2], "cint": batch_values[3],
        "wall_sec": batch_values[4], "timing": batch_values[5],
    }
    fixed = {
        "nw": fixed_values[0], "wave": fixed_values[1],
        "sint": fixed_values[2], "cint": fixed_values[3],
        "wall_sec": fixed_values[4], "timing": fixed_values[5],
    }
    flux, regular_wave, legacy_flux, batch_flux = flux_metrics(
        legacy, batch, 6561.0, 6564.2
    )
    result = {
        "model": "moderate_giant",
        "atmosphere_path": str(MODELS["moderate_giant"]),
        "window_A": [6561.0, 6564.2],
        "departure_grid": "nlte_H_pysme.grd",
        "input_lines": len(sme.linelist),
        "supported_lines": int(np.sum(supported)),
        "active_lines": int(np.sum(strong)),
        "legacy": {
            "points": legacy["nw"],
            "wall_sec": legacy["wall_sec"],
            "timing": legacy["timing"],
        },
        "batch": {
            "points": batch["nw"],
            "wall_sec": batch["wall_sec"],
            "timing": batch["timing"],
            "speedup": float(
                legacy["timing"].get("total_sec", legacy["wall_sec"])
                / batch["timing"].get("total_sec", batch["wall_sec"])
            ),
        },
        "grid_historical": grid_metrics(legacy["wave"], batch["wave"]),
        "historical_common_nodes": common_node_errors(legacy, batch),
        "controlled_fixed_reference": {
            "wave_bitwise_equal": bool(np.array_equal(fixed["wave"], batch["wave"])),
            "sint_bitwise_equal": bool(np.array_equal(fixed["sint"], batch["sint"])),
            "cint_bitwise_equal": bool(np.array_equal(fixed["cint"], batch["cint"])),
            "sint_max_abs": float(np.max(np.abs(fixed["sint"] - batch["sint"]))),
            "cint_max_abs": float(np.max(np.abs(fixed["cint"] - batch["cint"]))),
        },
        "flux_batch_minus_historical_legacy": flux,
    }
    arrays = {
        "nlte_halpha_legacy_wave": legacy["wave"],
        "nlte_halpha_legacy_sint": legacy["sint"],
        "nlte_halpha_legacy_cint": legacy["cint"],
        "nlte_halpha_batch_wave": batch["wave"],
        "nlte_halpha_batch_sint": batch["sint"],
        "nlte_halpha_batch_cint": batch["cint"],
        "nlte_halpha_regular_wave": regular_wave,
        "nlte_halpha_legacy_flux": legacy_flux,
        "nlte_halpha_batch_flux": batch_flux,
    }
    dll.SetLineInfoMode(0)
    dll.SetAdaptiveTransferGridMode("batched")
    return result, arrays


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--cases", nargs="+", choices=tuple(CASES), default=list(CASES))
    parser.add_argument(
        "--output", type=Path, default=ROOT / "analysis/spherical_batched_rkints_benchmark.json"
    )
    parser.add_argument(
        "--npz", type=Path, default=ROOT / "analysis/spherical_batched_rkints_benchmark.npz"
    )
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--nlte-only", action="store_true")
    args = parser.parse_args()

    os.environ["SME_TIMING"] = "1"
    if not DENSE_LINES.exists() and not args.nlte_only:
        raise FileNotFoundError(DENSE_LINES)
    linelists = {} if args.nlte_only else {"dense": ValdFile(str(DENSE_LINES))}
    if not args.nlte_only and any(CASES[name][3] == "halpha" for name in args.cases):
        linelists["halpha"] = ValdFile(str(ensure_halpha_lines()))

    output = {
        "setup": {
            "accrt": ACCRT,
            "physical_range_floor": RANGE_FLOOR,
            "accwi": ACCWI,
            "minimum_velocity_spacing_kms": DVEL_MIN_KMS,
            "line_lists": {
                "dense": str(DENSE_LINES),
                **(
                    {"halpha": str(ensure_halpha_lines())}
                    if "halpha" in linelists
                    else {}
                ),
            },
            "models": {key: str(value) for key, value in MODELS.items()},
            "nlte_empirically_validated": False,
        },
        "cases": {},
    }
    arrays = {}
    if args.resume and args.output.exists():
        output = json.loads(args.output.read_text())
    if args.resume and args.npz.exists():
        with np.load(args.npz) as old:
            arrays.update({key: old[key] for key in old.files})

    if args.nlte_only:
        result, nlte_arrays = run_spherical_halpha_nlte()
        output["spherical_halpha_nlte"] = result
        output["setup"]["nlte_empirically_validated"] = True
        arrays.update(nlte_arrays)
        args.output.write_text(json.dumps(json_clean(output), indent=2, sort_keys=True) + "\n")
        np.savez_compressed(args.npz, **arrays)
        print(
            f"spherical_halpha_nlte: legacy={result['legacy']['wall_sec']:.3f}s "
            f"batch={result['batch']['wall_sec']:.3f}s "
            f"speedup={result['batch']['speedup']:.2f}x",
            flush=True,
        )
        return

    for case in args.cases:
        print(f"running {case}", flush=True)
        result, case_arrays = run_case(case, linelists)
        output["cases"][case] = result
        arrays.update(case_arrays)
        args.output.write_text(json.dumps(json_clean(output), indent=2, sort_keys=True) + "\n")
        np.savez_compressed(args.npz, **arrays)
        print(
            f"{case}: legacy={result['legacy']['wall_sec']:.3f}s "
            f"batch={result['batch']['wall_sec']:.3f}s "
            f"speedup={result['batch']['speedup']:.2f}x",
            flush=True,
        )


if __name__ == "__main__":
    main()
