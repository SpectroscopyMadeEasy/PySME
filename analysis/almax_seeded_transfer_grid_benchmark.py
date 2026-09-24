"""Prototype and benchmark an ALMAX-seeded bounded transfer grid.

This is intentionally a diagnostic program.  It does not change the default
PySME/SMElib transfer-grid path.  Exact opacity/source and fine-grid transfer
solutions are evaluated once; the parameter sweep then operates on those
stored reference arrays before selected candidates are verified with a real
fixed-grid SMElib ``Transf`` call.
"""

from __future__ import annotations

import argparse
import copy
import json
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.ndimage import gaussian_filter1d

from pysme.abund import Abund
from pysme.linelist.vald import ValdFile
from pysme.sme import SME_Structure
from pysme.synthesize import Synthesizer


C_KMS = 299792.458
LO, HI = 5195.0, 5205.0
OUTPUT_STEP = 0.05
REFERENCE_STEP = 0.0003125
LINE_LIST = Path("/tmp/pysme_line_centric_windows/window_0010A_pad_150A.lin")
MODELS = {
    "solar_dwarf": (5777.0, 4.44, 0.0, 1.0),
    "metal_poor_dwarf": (6000.0, 4.0, -2.0, 1.2),
    "cool_metal_rich_dwarf": (4250.0, 4.5, 0.3, 1.0),
    "giant": (4500.0, 2.0, 0.0, 1.5),
}
ALMAX_EDGES = np.array([0.0, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, np.inf])
ALMAX_LABELS = ("<1e-6", "1e-6-1e-5", "1e-5-1e-4", "1e-4-1e-3", "1e-3-1e-2", ">=1e-2")


def make_sme(model_name: str, linelist: ValdFile, transfer_wave: np.ndarray) -> SME_Structure:
    teff, logg, monh, vmic = MODELS[model_name]
    output_wave = np.arange(LO, HI + 0.5 * OUTPUT_STEP, OUTPUT_STEP)
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
    sme.wave = [output_wave]
    sme.wint = [np.asarray(transfer_wave, dtype=float)]
    sme.vrad_flag = "none"
    sme.cscale_flag = "none"
    sme.accrt = 1e-4
    sme.accwi = 3e-3
    return sme


def integrate_mu_areas(mu, intensities):
    """Disk-integrate node intensities without wavelength resampling.

    This is the appropriate intrinsic/common-node diagnostic.  The production
    broadening integrator assumes a regular velocity grid and therefore must
    not be applied directly to an irregular adaptive transfer grid.
    """
    mu = np.asarray(mu, dtype=float)
    intensities = np.asarray(intensities, dtype=float)
    order = np.argsort(np.sqrt(1.0 - mu * mu))
    mu = mu[order]
    intensities = intensities[order]
    radius = np.sqrt(1.0 - mu * mu)
    if mu.size > 1:
        boundaries = np.sqrt(0.5 * (radius[:-1] ** 2 + radius[1:] ** 2))
        boundaries = np.concatenate(([0.0], boundaries, [1.0]))
        weights = boundaries[1:] ** 2 - boundaries[:-1] ** 2
    else:
        weights = np.ones(1, dtype=float)
    return np.pi * np.sum(weights[:, None] * intensities, axis=0)


def integrate_intrinsic(synth: Synthesizer, sme: SME_Structure, sint, cint):
    del synth
    line = integrate_mu_areas(sme.mu, sint)
    cont = integrate_mu_areas(sme.mu, cint)
    return np.asarray(line / cont), np.asarray(cont)


def direct_transf(synth: Synthesizer, sme: SME_Structure, wave: np.ndarray | None):
    dll = synth.get_dll()
    started = time.perf_counter()
    _, got_wave, sint, cint = dll.Transf(
        sme.mu,
        wave=wave,
        nwmax=40000 if wave is None else len(wave),
        accrt=sme.accrt,
        accwi=sme.accwi,
        keep_lineop=False,
        long_continuum=True,
    )
    elapsed = time.perf_counter() - started
    flux, continuum = integrate_intrinsic(synth, sme, sint, cint)
    return np.asarray(got_wave), flux, continuum, elapsed


def post_convolve(flux: np.ndarray, step: float, resolving_power: float) -> np.ndarray:
    fwhm = 0.5 * (LO + HI) / resolving_power
    sigma_pixels = fwhm / 2.354820045 / step
    return gaussian_filter1d(flux, sigma_pixels, mode="nearest")


def error_stats(candidate: np.ndarray, reference: np.ndarray) -> dict:
    delta = np.asarray(candidate) - np.asarray(reference)
    absolute = np.abs(delta)
    return {
        "max_abs": float(np.max(absolute)),
        "rms": float(np.sqrt(np.mean(delta * delta))),
        "p99_abs": float(np.quantile(absolute, 0.99)),
        "mean": float(np.mean(delta)),
    }


def is_molecular(species: np.ndarray) -> np.ndarray:
    elements = {
        "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg",
        "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr",
        "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr",
        "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
        "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd",
        "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf",
        "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po",
        "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U",
    }
    labels = np.asarray(species).astype(str)
    return np.array([label.split()[0] not in elements for label in labels], dtype=bool)


def legacy_seed_audit(wave, almax, strong, molecule, lo, hi):
    order = np.argsort(wave)
    wave = wave[order]
    almax = almax[order]
    strong = strong[order]
    molecule = molecule[order]
    selected = []
    last = lo
    for index in np.flatnonzero((wave > lo) & (wave < hi) & strong):
        min_spacing = wave[index] * 0.3 / C_KMS
        if wave[index] - last > min_spacing:
            selected.append(index)
            last = wave[index]
    selected = np.asarray(selected, dtype=int)
    bins = {}
    for label, left, right in zip(ALMAX_LABELS, ALMAX_EDGES[:-1], ALMAX_EDGES[1:]):
        mask = (almax[selected] >= left) & (almax[selected] < right)
        bins[label] = {
            "lines": int(np.sum(mask)),
            "initial_grid_points": int(2 * np.sum(mask)),
            "atomic_lines": int(np.sum(mask & ~molecule[selected])),
            "molecular_lines": int(np.sum(mask & molecule[selected])),
        }
    return {
        "accepted_line_centres": int(selected.size),
        "initial_line_driven_points": int(2 * selected.size),
        "atomic_centres": int(np.sum(~molecule[selected])),
        "molecular_centres": int(np.sum(molecule[selected])),
        "almax_bins": bins,
    }


def reference_index(wavelength: np.ndarray, fine_wave: np.ndarray) -> np.ndarray:
    index = np.rint((np.asarray(wavelength) - fine_wave[0]) / (fine_wave[1] - fine_wave[0])).astype(int)
    return np.clip(index, 0, fine_wave.size - 1)


def stage1_nodes(
    fine_wave,
    line_wave,
    almax,
    base_spacing,
    local_window,
    budget,
    position,
    strong_threshold,
):
    fine_step = fine_wave[1] - fine_wave[0]
    base_stride = max(1, int(round(base_spacing / fine_step)))
    base = np.arange(0, fine_wave.size, base_stride, dtype=int)
    if base[-1] != fine_wave.size - 1:
        base = np.append(base, fine_wave.size - 1)

    valid = (
        (line_wave >= fine_wave[0])
        & (line_wave <= fine_wave[-1])
        & np.isfinite(almax)
        & (almax > 0)
    )
    lw = line_wave[valid]
    la = almax[valid]
    order = np.argsort(lw)
    lw, la = lw[order], la[order]

    seed = []
    start = fine_wave[0]
    while start < fine_wave[-1]:
        stop = min(start + local_window, fine_wave[-1])
        left = np.searchsorted(lw, start, side="left")
        right = np.searchsorted(lw, stop, side="left")
        if right > left and float(np.sum(la[left:right])) >= budget:
            if position == "centroid":
                loc = float(np.average(lw[left:right], weights=la[left:right]))
            else:
                loc = 0.5 * (start + stop)
            seed.append(int(reference_index(np.array([loc]), fine_wave)[0]))
        start = stop

    strong = reference_index(lw[la >= strong_threshold], fine_wave)
    nodes = np.unique(np.concatenate((base, np.asarray(seed, dtype=int), strong)))
    return nodes, {
        "base_points": int(base.size),
        "importance_seed_points": int(np.setdiff1d(seed, base).size),
        "strong_line_seed_points": int(np.setdiff1d(strong, np.union1d(base, seed)).size),
        "local_bins_above_budget": int(len(seed)),
    }


def quantity_error(exact, predicted, floor_fraction=1e-8):
    scale = max(float(np.max(np.abs(exact))) * floor_fraction, np.finfo(float).tiny)
    return float(np.max(np.abs(exact - predicted) / np.maximum(np.abs(exact), scale)))


@dataclass
class RefinementResult:
    nodes: np.ndarray
    added: int
    unresolved: int
    worst_unresolved: float
    min_spacing: float
    max_depth_reached: int


def refine_nodes(
    initial_nodes,
    chi,
    eta,
    line_indices,
    line_almax,
    tolerance,
    probe_strategy,
    probe_threshold,
    min_steps=1,
    max_points=20000,
    max_depth=16,
):
    line_order = np.argsort(line_indices)
    line_indices = np.asarray(line_indices, dtype=int)[line_order]
    line_almax = np.asarray(line_almax, dtype=float)[line_order]
    accepted = set(int(value) for value in initial_nodes)
    stack = [(int(a), int(b), 0) for a, b in zip(initial_nodes[:-1], initial_nodes[1:])]
    unresolved = 0
    worst_unresolved = 0.0
    max_depth_reached = 0

    while stack:
        left, right, depth = stack.pop()
        if right - left <= 1:
            continue
        max_depth_reached = max(max_depth_reached, depth)
        fraction_denom = float(right - left)
        probes = [left + (right - left) // 2]
        ll = np.searchsorted(line_indices, left + 1, side="left")
        rr = np.searchsorted(line_indices, right, side="left")
        if rr > ll and probe_strategy != "midpoint":
            local_indices = line_indices[ll:rr]
            local_almax = line_almax[ll:rr]
            if probe_strategy == "strongest":
                best = int(np.argmax(local_almax))
                if local_almax[best] >= probe_threshold:
                    probes.append(int(local_indices[best]))
            elif probe_strategy == "all_above":
                probes.extend(int(value) for value in local_indices[local_almax >= probe_threshold])
        probes = sorted(set(value for value in probes if left < value < right))

        worst_error = -1.0
        worst_probe = probes[0]
        for probe in probes:
            fraction = (probe - left) / fraction_denom
            chi_interp = chi[left] + fraction * (chi[right] - chi[left])
            eta_interp = eta[left] + fraction * (eta[right] - eta[left])
            error = max(
                quantity_error(chi[probe], chi_interp),
                quantity_error(eta[probe], eta_interp),
            )
            if error > worst_error:
                worst_error = error
                worst_probe = probe

        if worst_error <= tolerance:
            continue
        if (
            right - left <= min_steps
            or depth >= max_depth
            or len(accepted) >= max_points
        ):
            unresolved += 1
            worst_unresolved = max(worst_unresolved, worst_error)
            continue
        accepted.add(worst_probe)
        stack.append((left, worst_probe, depth + 1))
        stack.append((worst_probe, right, depth + 1))

    nodes = np.asarray(sorted(accepted), dtype=int)
    return RefinementResult(
        nodes=nodes,
        added=int(nodes.size - len(initial_nodes)),
        unresolved=unresolved,
        worst_unresolved=worst_unresolved,
        min_spacing=float(np.min(np.diff(nodes))),
        max_depth_reached=max_depth_reached,
    )


def feature_windows(line_wave, almax, molecule):
    inside = (line_wave >= LO) & (line_wave <= HI) & np.isfinite(almax)
    atomic = inside & ~molecule
    molecular = inside & molecule
    weak_candidates = np.flatnonzero(atomic & (almax >= 1e-4) & (almax < 1e-2))
    if weak_candidates.size:
        weak = weak_candidates[np.argmin(np.abs(np.log10(almax[weak_candidates]) + 3.0))]
    else:
        weak = np.flatnonzero(atomic)[np.argmin(almax[atomic])]
    strong = np.flatnonzero(atomic)[np.argmax(almax[atomic])]
    mol_wave = line_wave[molecular]
    mol_almax = almax[molecular]
    centres = np.arange(LO + 0.1, HI - 0.1, 0.02)
    sums = np.array([np.sum(mol_almax[np.abs(mol_wave - centre) <= 0.1]) for centre in centres])
    forest = float(centres[np.argmax(sums)])
    return {
        "weak_atomic": (float(line_wave[weak] - 0.08), float(line_wave[weak] + 0.08)),
        "strong_atomic": (float(line_wave[strong] - 0.15), float(line_wave[strong] + 0.15)),
        "molecular_forest": (forest - 0.15, forest + 0.15),
        "worst_coarse_feature": (5200.75, 5201.05),
    }


def equivalent_widths(wave, flux, windows):
    result = {}
    for name, (left, right) in windows.items():
        mask = (wave >= left) & (wave <= right)
        result[name] = float(np.trapezoid(1.0 - flux[mask], wave[mask]))
    return result


def evaluate_candidate(nodes, fine_wave, reference_flux, windows):
    wave = fine_wave[nodes]
    interpolated = np.interp(fine_wave, wave, reference_flux[nodes])
    result = {
        "transfer_points": int(nodes.size),
        "points_per_A": float(nodes.size / (fine_wave[-1] - fine_wave[0])),
        "minimum_spacing_A": float(np.min(np.diff(wave))),
        "intrinsic": error_stats(interpolated, reference_flux),
        "ew": equivalent_widths(fine_wave, interpolated, windows),
        "post_convolution": {},
    }
    for resolving_power in (20000.0, 60000.0):
        ref_conv = post_convolve(reference_flux, fine_wave[1] - fine_wave[0], resolving_power)
        got_conv = post_convolve(interpolated, fine_wave[1] - fine_wave[0], resolving_power)
        result["post_convolution"][str(int(resolving_power))] = error_stats(got_conv, ref_conv)
    return result, interpolated


def run_model(model_name, linelist, include_legacy=False):
    coarse = np.arange(LO, HI + 0.5 * OUTPUT_STEP, OUTPUT_STEP)
    sme = make_sme(model_name, linelist, coarse)
    synth = Synthesizer()
    started = time.perf_counter()
    synth.synthesize_spectrum(
        sme,
        updateStructure=False,
        reuse_wavelength_grid=True,
        radial_velocity_mode="fast",
        linelist_mode="all",
        smelib_lineinfo_mode=0,
    )
    setup_time = time.perf_counter() - started
    dll = synth.get_dll()

    line_wave = np.asarray(sme.linelist["wlcent"], dtype=float)
    almax = np.asarray(sme.linelist["almax_ratio"], dtype=float)
    strong = np.asarray(sme.linelist["strong"], dtype=bool)
    molecule = is_molecular(np.asarray(sme.linelist["species"]))

    legacy = None
    legacy_arrays = {}
    if include_legacy:
        legacy_wave, legacy_flux, legacy_cont, legacy_sec = direct_transf(synth, sme, None)
        audit = legacy_seed_audit(
            line_wave,
            almax,
            strong,
            molecule,
            float(legacy_wave[0]),
            float(legacy_wave[-1]),
        )
        audit["refinement_or_other_points"] = int(
            legacy_wave.size - audit["initial_line_driven_points"] - 2
        )
        legacy = {
            "raw_transfer_points": int(legacy_wave.size),
            "transf_sec": legacy_sec,
            "wave_min": float(legacy_wave[0]),
            "wave_max": float(legacy_wave[-1]),
            "minimum_spacing_A": float(np.min(np.diff(legacy_wave))),
            "seed_audit": audit,
        }
        legacy_arrays = {
            f"{model_name}_legacy_wave": legacy_wave,
            f"{model_name}_legacy_flux": legacy_flux,
            f"{model_name}_legacy_continuum": legacy_cont,
        }

    fine_wave = np.arange(LO, HI + 0.5 * REFERENCE_STEP, REFERENCE_STEP)
    ref_wave, ref_flux, ref_cont, ref_sec = direct_transf(synth, sme, fine_wave)
    opacity_started = time.perf_counter()
    chi, continuum_chi, scatter, source, continuum_source = dll.GetOpacityAtWaves(fine_wave)
    opacity_sec = time.perf_counter() - opacity_started
    eta = chi * source
    _, repeat_flux, repeat_cont, repeat_sec = direct_transf(synth, sme, fine_wave)
    reference_repeat_error = error_stats(repeat_flux, ref_flux)
    ref_flux = repeat_flux
    ref_cont = repeat_cont

    valid_lines = (
        (line_wave >= LO)
        & (line_wave <= HI)
        & np.isfinite(almax)
        & strong
    )
    line_indices = reference_index(line_wave[valid_lines], fine_wave)
    line_almax = almax[valid_lines]
    windows = feature_windows(line_wave, almax, molecule)
    reference_ew = equivalent_widths(fine_wave, ref_flux, windows)

    result = {
        "model": model_name,
        "setup_sec": setup_time,
        "reference": {
            "spacing_A": REFERENCE_STEP,
            "points": int(fine_wave.size),
            "transf_sec": ref_sec,
            "repeat_transf_sec": repeat_sec,
            "repeat_after_opacity_error": reference_repeat_error,
            "opacity_source_sec": opacity_sec,
            "ew": reference_ew,
        },
        "line_counts": {
            "inside": int(np.sum((line_wave >= LO) & (line_wave <= HI))),
            "selected_inside": int(np.sum(valid_lines)),
            "atomic_inside": int(np.sum(valid_lines & ~molecule)),
            "molecular_inside": int(np.sum(valid_lines & molecule)),
        },
        "legacy": legacy,
        "fixed": {},
        "stage1_sweep": [],
        "focused_sweep": [],
        "verified_candidates": [],
        "feature_windows": {name: list(bounds) for name, bounds in windows.items()},
    }

    for spacing in (0.05, 0.02, 0.01, 0.005, 0.0025, 0.00125, 0.000625):
        stride = max(1, int(round(spacing / REFERENCE_STEP)))
        nodes = np.arange(0, fine_wave.size, stride, dtype=int)
        if nodes[-1] != fine_wave.size - 1:
            nodes = np.append(nodes, fine_wave.size - 1)
        record, _ = evaluate_candidate(nodes, fine_wave, ref_flux, windows)
        result["fixed"][str(spacing)] = record

    sweep_configs = []
    for base_spacing in (0.05, 0.02, 0.01):
        for local_window in (0.02, 0.05, 0.10):
            for budget in (1e-2, 3e-3, 1e-3, 3e-4):
                sweep_configs.append((base_spacing, local_window, budget, "centre", 10.0, 1e-3, "strongest", 1e-3))

    def one_config(config):
        base_spacing, local_window, budget, position, strong_threshold, tolerance, strategy, probe_threshold = config
        initial, counts = stage1_nodes(
            fine_wave,
            line_wave,
            almax,
            base_spacing,
            local_window,
            budget,
            position,
            strong_threshold,
        )
        refined = refine_nodes(
            initial,
            chi,
            eta,
            line_indices,
            line_almax,
            tolerance,
            strategy,
            probe_threshold,
        )
        metrics, interpolated = evaluate_candidate(refined.nodes, fine_wave, ref_flux, windows)
        ew_delta = {name: metrics["ew"][name] - reference_ew[name] for name in windows}
        return {
            "config": {
                "base_spacing_A": base_spacing,
                "local_window_A": local_window,
                "importance_budget": budget,
                "seed_position": position,
                "strong_line_threshold": strong_threshold,
                "interpolation_tolerance": tolerance,
                "probe_strategy": strategy,
                "probe_threshold": probe_threshold,
            },
            "counts": {
                **counts,
                "error_refinement_points": refined.added,
                "final_transfer_points": int(refined.nodes.size),
                "unresolved_intervals": refined.unresolved,
                "worst_unresolved_error": refined.worst_unresolved,
                "max_refinement_depth": refined.max_depth_reached,
            },
            "metrics": metrics,
            "ew_delta_A": ew_delta,
            "_nodes": refined.nodes,
            "_interpolated": interpolated,
        }

    for config in sweep_configs:
        record = one_config(config)
        result["stage1_sweep"].append({key: value for key, value in record.items() if not key.startswith("_")})

    ranked = sorted(
        zip(sweep_configs, result["stage1_sweep"]),
        key=lambda item: (
            item[1]["counts"]["unresolved_intervals"] > 0,
            item[1]["metrics"]["intrinsic"]["max_abs"] > 1e-3,
            item[1]["counts"]["final_transfer_points"],
            item[1]["metrics"]["intrinsic"]["max_abs"],
        ),
    )
    best = ranked[0][0]
    base_spacing, local_window, budget, _, _, _, _, _ = best
    focused = []
    for position in ("centre", "centroid"):
        for strong_threshold in (1.0, 10.0, 100.0):
            for tolerance in (3e-3, 1e-3, 3e-4, 1e-4):
                for strategy, probe_threshold in (("midpoint", 0.0), ("strongest", 1e-3), ("all_above", 10.0)):
                    focused.append((base_spacing, local_window, budget, position, strong_threshold, tolerance, strategy, probe_threshold))

    focused_full = []
    for config in focused:
        record = one_config(config)
        focused_full.append(record)
        result["focused_sweep"].append({key: value for key, value in record.items() if not key.startswith("_")})

    verify = sorted(
        focused_full,
        key=lambda item: (
            item["counts"]["unresolved_intervals"] > 0,
            item["metrics"]["intrinsic"]["max_abs"] > 3e-4,
            item["counts"]["final_transfer_points"],
        ),
    )[:2]
    candidate_arrays = {}
    for number, candidate in enumerate(verify):
        nodes = candidate["_nodes"]
        wave = fine_wave[nodes]
        got_wave, got_flux_nodes, got_cont_nodes, elapsed = direct_transf(synth, sme, wave)
        got_flux = np.interp(fine_wave, got_wave, got_flux_nodes)
        verification = {
            "candidate": number,
            "config": candidate["config"],
            "transfer_points": int(nodes.size),
            "transf_sec": elapsed,
            "common_node_intrinsic": error_stats(
                got_flux_nodes, ref_flux[nodes]
            ),
            "actual_intrinsic": error_stats(got_flux, ref_flux),
            "simulation_vs_actual": error_stats(candidate["_interpolated"], got_flux),
        }
        result["verified_candidates"].append(verification)
        candidate_arrays[f"{model_name}_candidate_{number}_wave"] = got_wave
        candidate_arrays[f"{model_name}_candidate_{number}_flux"] = got_flux_nodes

    arrays = {
        f"{model_name}_reference_wave": ref_wave,
        f"{model_name}_reference_flux": ref_flux,
        f"{model_name}_reference_continuum": ref_cont,
        f"{model_name}_chi_total": chi,
        f"{model_name}_eta_total": eta,
        f"{model_name}_continuum_chi": continuum_chi,
        f"{model_name}_scatter": scatter,
        f"{model_name}_continuum_source": continuum_source,
        f"{model_name}_line_wave": line_wave[valid_lines],
        f"{model_name}_line_almax": almax[valid_lines],
        f"{model_name}_line_molecular": molecule[valid_lines].astype(np.uint8),
        **legacy_arrays,
        **candidate_arrays,
    }
    return result, arrays


def json_clean(value):
    if isinstance(value, dict):
        return {str(key): json_clean(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_clean(item) for item in value]
    if isinstance(value, np.generic):
        return value.item()
    return value


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--models", nargs="+", choices=tuple(MODELS), default=list(MODELS))
    parser.add_argument("--legacy", action="store_true")
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--output", type=Path, default=Path("analysis/almax_seeded_transfer_grid_benchmark.json"))
    parser.add_argument("--npz", type=Path, default=Path("analysis/almax_seeded_transfer_grid_benchmark.npz"))
    args = parser.parse_args()
    if not LINE_LIST.exists():
        raise FileNotFoundError(LINE_LIST)
    linelist = ValdFile(str(LINE_LIST))
    if args.resume and args.output.exists():
        output = json.loads(args.output.read_text())
        output["setup"]["models"] = sorted(set(output["setup"].get("models", [])) | set(args.models))
    else:
        output = {
            "setup": {
                "window_A": [LO, HI],
                "output_step_A": OUTPUT_STEP,
                "reference_step_A": REFERENCE_STEP,
                "line_list": str(LINE_LIST),
                "input_lines": int(len(linelist)),
                "models": args.models,
            },
            "models": {},
        }
    arrays = {}
    if args.resume and args.npz.exists():
        with np.load(args.npz) as existing:
            arrays.update({name: existing[name] for name in existing.files})
    for model_name in args.models:
        print(f"running {model_name}", flush=True)
        result, model_arrays = run_model(
            model_name,
            linelist,
            include_legacy=args.legacy and model_name == "cool_metal_rich_dwarf",
        )
        output["models"][model_name] = result
        arrays.update(model_arrays)
        args.output.write_text(json.dumps(json_clean(output), indent=2, sort_keys=True) + "\n")
        np.savez_compressed(args.npz, **arrays)
    print(args.output)
    print(args.npz)


if __name__ == "__main__":
    main()
