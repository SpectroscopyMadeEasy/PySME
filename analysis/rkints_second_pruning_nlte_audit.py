"""Append a small H-alpha NLTE sanity check to the RKINTS pruning audit."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from pysme.synthesize import Synthesizer

from rkints_second_pruning_audit import (
    FINE_STEP,
    LOW_RANGE_FLOOR,
    ROOT,
    compact_indices,
    compute_low_floor_lineinfo,
    input_line_state,
    json_clean,
    metrics_on_common_grid,
    normalized_intrinsic,
    raw_transf,
    read_prune_log,
    run_prune_diagnostic,
    timed_fixed_transf,
)


def main():
    # Reuse the tested regression fixture so the departure coefficients and
    # transition matching are identical to PySME's H-alpha NLTE regression.
    import sys

    sys.path.insert(0, str(ROOT))
    from test.test_halpha_regression import _make_structure

    sme = _make_structure("sun_nlte")
    sme.vmac = 0.0
    sme.vsini = 0.0
    sme.ipres = 1e7
    sme.accrt = 1e-4
    sme.accwi = 3e-3
    synth = Synthesizer()
    synth.synthesize_spectrum(sme)
    dll = synth.get_dll()

    log_path = ROOT / "analysis" / "rkints_second_pruning_logs" / "sun_halpha_nlte.csv"
    legacy_wave, legacy_flux, legacy_sec, rows = run_prune_diagnostic(
        dll, sme, "immediate", log_path
    )
    seed_rows = [row for row in rows if row["phase"] == "seed"]
    removed_rows = [row for row in seed_rows if row["applied"]]
    removed_indices = np.asarray([row["line_index"] for row in removed_rows], dtype=int)
    nlte_flags = np.asarray(sme.nlte.flags, dtype=bool)

    prod_mask = np.asarray(sme.linelist["strong"], dtype=bool)
    prod_ranges = np.column_stack(
        (
            np.asarray(sme.linelist["line_range_s"], dtype=float),
            np.asarray(sme.linelist["line_range_e"], dtype=float),
        )
    )
    final_mask = prod_mask.copy()
    final_mask[removed_indices] = False
    wfirst, wlast = float(legacy_wave[0]), float(legacy_wave[-1])
    fine_wave = np.arange(wfirst, wlast + 0.5 * FINE_STEP, FINE_STEP)
    low_almax, low_ranges = compute_low_floor_lineinfo(dll, sme, wfirst, wlast)
    supported = ~np.asarray(sme.line_ion_mask, dtype=bool)
    valid_low = supported & np.isfinite(low_almax)
    modes = {
        "L_final_mask": (prod_ranges, final_mask),
        "A_almax_only": (prod_ranges, prod_mask),
        "P_1e-5": (low_ranges, valid_low & (low_almax >= 1e-5)),
        "F_1e-6": (low_ranges, valid_low & (low_almax >= LOW_RANGE_FLOOR)),
    }
    fluxes = {}
    timings = {}
    for name, (ranges, mask) in modes.items():
        _, flux, timing = timed_fixed_transf(dll, sme, fine_wave, ranges, mask)
        fluxes[name] = flux
        timings[name] = timing

    lo, hi = 6561.0, 6564.2
    centre = (fine_wave >= lo) & (fine_wave <= hi)
    wave = fine_wave[centre]
    reference = fluxes["F_1e-6"][centre]
    windows = {"halpha": (lo, hi)}
    comparisons = {
        name: metrics_on_common_grid(flux[centre], reference, wave, windows)
        for name, flux in fluxes.items()
    }
    comparisons["L_historical"] = metrics_on_common_grid(
        np.interp(wave, legacy_wave, legacy_flux), reference, wave, windows
    )

    result = {
        "case": "sun_halpha_nlte",
        "window_A": [lo, hi],
        "input_lines": int(len(sme.linelist)),
        "supported_lines": int(compact_indices(sme).size),
        "almax_lines": int(np.sum(prod_mask)),
        "nlte_lines": int(np.sum(nlte_flags)),
        "nlte_line_indices": np.flatnonzero(nlte_flags),
        "seed_centres_tested": int(len(seed_rows)),
        "additionally_mark2": int(len(removed_rows)),
        "nlte_lines_marked2": int(
            np.sum(nlte_flags[removed_indices]) if removed_indices.size else 0
        ),
        "legacy_transf_sec": legacy_sec,
        "mode_line_counts": {
            name: int(np.sum(mask)) for name, (_, mask) in modes.items()
        },
        "comparisons_to_F_1e-6": comparisons,
        "fine_grid_timings": timings,
        "source_semantics": (
            "OPMTRX checks MARK before its NLTE branch; MARK=2 therefore skips "
            "both NLTE extinction and the line emissivity/source contribution"
        ),
    }

    json_path = ROOT / "analysis" / "rkints_second_pruning_audit.json"
    npz_path = ROOT / "analysis" / "rkints_second_pruning_audit.npz"
    output = json.loads(json_path.read_text())
    output["nlte_sanity"] = result
    json_path.write_text(json.dumps(json_clean(output), indent=2, sort_keys=True) + "\n")

    arrays = {}
    with np.load(npz_path) as existing:
        arrays.update({name: existing[name] for name in existing.files})
    arrays["nlte_sun_halpha_legacy_wave"] = legacy_wave
    arrays["nlte_sun_halpha_legacy_flux"] = legacy_flux
    arrays["nlte_sun_halpha_fine_wave"] = fine_wave
    for name, flux in fluxes.items():
        arrays[f"nlte_sun_halpha_{name}_flux"] = flux
    np.savez_compressed(npz_path, **arrays)
    print(json.dumps(json_clean(result), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
