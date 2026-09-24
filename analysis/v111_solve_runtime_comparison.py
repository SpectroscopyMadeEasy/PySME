"""End-to-end abundance ``solve`` benchmark across PySME revisions.

Run this script once in each isolated installation.  The target spectrum is a
precomputed, noise-free spectrum made with the same revision and its defaults;
target construction is therefore excluded from the timed solve.  Starting
0.10 dex away from the target keeps the optimization problem self-consistent
in both revisions, while the synthesis-call counter makes different optimizer
iteration counts visible.
"""

from __future__ import annotations

import argparse
import copy
import json
import platform
import time
from pathlib import Path

import numpy as np

import pysme
from pysme.abund import Abund
from pysme.linelist.vald import ValdFile
from pysme.sme import MASK_VALUES, SME_Structure
from pysme.solve import solve
from pysme.synthesize import Synthesizer

from v111_runtime_comparison import CENTER_A, MODELS, STEP_A, window_path


def make_sme(
    model: tuple,
    linelist: ValdFile,
    width: float,
    transfer_grid: str,
) -> SME_Structure:
    _, teff, logg, monh, vmic = model
    lo = CENTER_A - width / 2
    hi = CENTER_A + width / 2
    wave = np.linspace(lo, hi, int(round(width / STEP_A)) + 1)

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
    sme.wave = [wave]
    if transfer_grid == "fixed":
        sme.wint = [wave]
    sme.vrad_flag = "none"
    sme.cscale_flag = "none"
    sme.accrt = 1e-4
    return sme


def build_reachable_target(
    args: argparse.Namespace,
    model: tuple,
    linelist: ValdFile,
) -> None:
    """Build a target after exercising the same abundance-change lifecycle."""
    sme = make_sme(model, linelist, args.width, args.transfer_grid)
    fe_index = sme.abund.elem_dict["Fe"]
    target_fe = float(sme.abund.get_pattern(type="H=12", raw=True)[fe_index])
    sme["Abund Fe"] = target_fe + args.initial_offset
    synth = Synthesizer()
    common = {
        "updateStructure": False,
        "reuse_wavelength_grid": args.transfer_grid == "fixed",
        "radial_velocity_mode": "fast",
        "linelist_mode": "all",
        "smelib_lineinfo_mode": 0,
    }
    synth.synthesize_spectrum(sme, passLineList=True, passAtmosphere=True, **common)
    sme["Abund Fe"] = target_fe
    result = synth.synthesize_spectrum(
        sme,
        passLineList=True,
        passAtmosphere=True,
        **common,
    )
    target = np.asarray(result[1][0], dtype=np.float64).copy()
    np.savez_compressed(args.build_target, target=target)
    print(
        json.dumps(
            {
                "version": str(pysme.__version__),
                "model": args.model,
                "output": str(args.build_target),
                "target_key": "target",
                "selection_recomputed_after_offset_dex": args.initial_offset,
                "transfer_grid": args.transfer_grid,
                "target_fe_pattern": target_fe,
                "wave_points": int(target.size),
                "input_lines": int(len(linelist)),
            },
            indent=2,
            sort_keys=True,
        ),
        flush=True,
    )


def run(args: argparse.Namespace) -> None:
    models = {model[0]: model for model in MODELS}
    model = models[args.model]
    linelist_path = window_path(args.width)

    parse_start = time.perf_counter()
    linelist = ValdFile(str(linelist_path))
    parse_sec = time.perf_counter() - parse_start
    if args.build_target is not None:
        build_reachable_target(args, model, linelist)
        return

    sme = make_sme(model, linelist, args.width, args.transfer_grid)

    with np.load(args.target) as target_data:
        target = np.asarray(target_data[args.target_key], dtype=np.float64).copy()
    if target.shape != np.asarray(sme.wave[0]).shape:
        raise ValueError(
            f"target shape {target.shape} does not match wave {sme.wave[0].shape}"
        )

    sme.spec = [target]
    sme.uncs = [np.full(target.size, args.uncertainty, dtype=np.float64)]
    sme.mask = [np.full(target.size, MASK_VALUES.LINE, dtype=np.int16)]

    fe_index = sme.abund.elem_dict["Fe"]
    target_fe_pattern = float(sme.abund.get_pattern(type="H=12", raw=True)[fe_index])
    sme["Abund Fe"] = target_fe_pattern + args.initial_offset

    synthesis_calls = 0
    synthesis_sec = 0.0
    original_synthesize = Synthesizer.synthesize_spectrum

    def counted_synthesize(self, *call_args, **call_kwargs):
        nonlocal synthesis_calls, synthesis_sec
        started = time.perf_counter()
        try:
            return original_synthesize(self, *call_args, **call_kwargs)
        finally:
            synthesis_calls += 1
            synthesis_sec += time.perf_counter() - started

    Synthesizer.synthesize_spectrum = counted_synthesize
    started = time.perf_counter()
    try:
        result = solve(
            sme,
            ["Abund Fe"],
            filename=None,
            restore=False,
            linelist_mode="all",
            smelib_lineinfo_mode=0,
        )
    finally:
        solve_sec = time.perf_counter() - started
        Synthesizer.synthesize_spectrum = original_synthesize

    fitted_pattern = float(result.fitresults.values[0])
    fitted_effective = fitted_pattern + float(result.monh)
    target_effective = target_fe_pattern + float(result.monh)
    final_flux = np.asarray(result.synth[0], dtype=np.float64)
    residual = final_flux - target

    record = {
        "version": str(pysme.__version__),
        "python": platform.python_version(),
        "platform": platform.platform(),
        "model": args.model,
        "width_A": args.width,
        "center_A": CENTER_A,
        "step_A": STEP_A,
        "transfer_grid": args.transfer_grid,
        "wave_points": int(target.size),
        "input_lines": int(len(linelist)),
        "accrt": float(sme.accrt),
        "target_file": str(args.target),
        "target_key": args.target_key,
        "target_kind": "same-revision default, noise-free, precomputed",
        "uncertainty": args.uncertainty,
        "initial_offset_dex": args.initial_offset,
        "target_fe_pattern": target_fe_pattern,
        "initial_fe_pattern": target_fe_pattern + args.initial_offset,
        "fitted_fe_pattern": fitted_pattern,
        "target_fe_effective": target_effective,
        "fitted_fe_effective": fitted_effective,
        "fit_error_dex": fitted_pattern - target_fe_pattern,
        "solve_sec": solve_sec,
        "synthesis_calls": synthesis_calls,
        "synthesis_sec": synthesis_sec,
        "non_synthesis_sec": solve_sec - synthesis_sec,
        "sec_per_synthesis_call": synthesis_sec / synthesis_calls,
        "residual_iterations": int(result.fitresults.iterations),
        "chisq": float(result.fitresults.chisq),
        "final_flux_max_abs_residual": float(np.max(np.abs(residual))),
        "final_flux_rms_residual": float(np.sqrt(np.mean(residual**2))),
        "linelist_parse_sec_excluded": parse_sec,
    }
    args.output.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
    print(json.dumps(record, indent=2, sort_keys=True), flush=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--model",
        choices=[model[0] for model in MODELS],
        required=True,
    )
    parser.add_argument("--width", type=float, default=10.0)
    parser.add_argument(
        "--transfer-grid",
        choices=("fixed", "wave-only"),
        default="fixed",
    )
    parser.add_argument("--target", type=Path)
    parser.add_argument("--target-key")
    parser.add_argument("--build-target", type=Path)
    parser.add_argument("--initial-offset", type=float, default=0.10)
    parser.add_argument("--uncertainty", type=float, default=0.01)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.build_target is None:
        if args.target is None or args.target_key is None or args.output is None:
            parser.error("--target, --target-key, and --output are required for solve")
    return args


if __name__ == "__main__":
    run(parse_args())
