"""Check the metal-poor full-reference floor from 1e-6 to 1e-7."""

from __future__ import annotations

import copy
import json

import numpy as np

from almax_seeded_transfer_grid_benchmark import LINE_LIST
from pysme.linelist.vald import ValdFile

from rkints_second_pruning_audit import (
    ROOT,
    compact_indices,
    compute_low_floor_lineinfo,
    error_stats,
    json_clean,
    setup_model,
    timed_fixed_transf,
)


def main():
    json_path = ROOT / "analysis" / "rkints_second_pruning_audit.json"
    npz_path = ROOT / "analysis" / "rkints_second_pruning_audit.npz"
    output = json.loads(json_path.read_text())
    with np.load(npz_path) as existing:
        arrays = {name: existing[name] for name in existing.files}
    fine_wave = arrays["metal_poor_dwarf_fine_wave"]
    flux_1e6 = arrays["metal_poor_dwarf_F_1e-6_fine_flux"]

    linelist = ValdFile(str(LINE_LIST))
    sme, synth, _ = setup_model("metal_poor_dwarf", copy.deepcopy(linelist))
    dll = synth.get_dll()
    almax_1e7, ranges_1e7 = compute_low_floor_lineinfo(
        dll, sme, float(fine_wave[0]), float(fine_wave[-1]), floor=1e-7
    )
    supported = ~np.asarray(sme.line_ion_mask, dtype=bool)
    mask_1e7 = supported & np.isfinite(almax_1e7) & (almax_1e7 >= 1e-7)
    _, flux_1e7, timing = timed_fixed_transf(
        dll, sme, fine_wave, ranges_1e7, mask_1e7
    )
    centre = (fine_wave >= 5195.0) & (fine_wave <= 5205.0)
    convergence = {
        "floor_1e-6_lines": int(
            output["models"]["metal_poor_dwarf"]["line_counts"]["full_1e-6"]
        ),
        "floor_1e-7_lines": int(np.sum(mask_1e7)),
        "F_1e-6_minus_F_1e-7": error_stats(
            flux_1e6[centre], flux_1e7[centre]
        ),
        "F_1e-7_timing": timing,
    }
    output["models"]["metal_poor_dwarf"]["reference_convergence"] = convergence
    json_path.write_text(json.dumps(json_clean(output), indent=2, sort_keys=True) + "\n")
    arrays["metal_poor_dwarf_F_1e-7_fine_flux"] = flux_1e7
    arrays["metal_poor_dwarf_almax_1e-7"] = almax_1e7
    np.savez_compressed(npz_path, **arrays)
    print(json.dumps(json_clean(convergence), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
