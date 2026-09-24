"""Version-to-version runtime benchmark for PySME v1.1.1 and current HEAD.

Run this script in an isolated installation of each revision.  It uses the
same extracted VALD files, atmospheres, fixed transfer grid, abundance trial
sequence, and machine for both versions.  VALD parsing is reported but excluded
from synthesis timings.

The default profile intentionally uses each revision's defaults: v1.1.1 uses
internal line selection and exact continuum opacity, while current HEAD uses
ALMAX selection and the adaptive continuum grid.  The ``reference`` profile
requests the legacy internal/exact path where the installed revision supports
the corresponding controls.
"""

from __future__ import annotations

import argparse
import copy
import inspect
import json
import time
from pathlib import Path

import numpy as np

import pysme
from pysme.abund import Abund
from pysme.linelist.vald import ValdFile
from pysme.sme import SME_Structure
from pysme.synthesize import Synthesizer


WINDOW_DIR = Path("/tmp/pysme_line_centric_windows")
STEP_A = 0.05
CENTER_A = 5200.0
FE_OFFSETS_DEX = (0.00, 0.05, -0.05, 0.10, -0.10, 0.00)
MODELS = (
    ("solar_dwarf", 5777.0, 4.44, 0.0, 1.0),
    ("metal_poor_dwarf", 6000.0, 4.0, -2.0, 1.2),
    ("cool_metal_rich_dwarf", 4250.0, 4.5, 0.3, 1.0),
)


def window_path(width: float) -> Path:
    return WINDOW_DIR / f"window_{int(width):04d}A_pad_150A.lin"


def make_sme(
    model: tuple,
    linelist: ValdFile,
    width: float,
    profile: str,
    transfer_grid: str,
):
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

    if profile == "reference":
        if hasattr(sme, "line_select_method"):
            sme.line_select_method = "internal"
        if hasattr(sme, "continuum_grid"):
            sme.continuum_grid = "exact"
    return sme


def synthesize(synth: Synthesizer, sme: SME_Structure, **kwargs):
    started = time.perf_counter()
    result = synth.synthesize_spectrum(
        sme,
        updateStructure=False,
        reuse_wavelength_grid=kwargs.pop("reuse_wavelength_grid", True),
        radial_velocity_mode="fast",
        linelist_mode="all",
        smelib_lineinfo_mode=0,
        **kwargs,
    )
    wall = time.perf_counter() - started
    flux = np.asarray(result[1][0], dtype=np.float64).copy()
    return wall, flux


def run_case(
    model: tuple,
    linelist: ValdFile,
    width: float,
    profile: str,
    offsets: tuple[float, ...],
    eos_warm: bool,
    transfer_grid: str,
):
    sme = make_sme(model, linelist, width, profile, transfer_grid)
    synth = Synthesizer()
    dll = synth.get_dll()
    supports_prepared = "passAbund" in inspect.signature(
        synth.synthesize_spectrum
    ).parameters
    warm_control = getattr(dll, "SetEosWarmStartMode", None)
    if warm_control is not None and eos_warm:
        warm_control(True)

    base_fe = float(sme.abund.A["Fe"])
    initial_sec, initial_flux = synthesize(
        synth,
        sme,
        passLineList=True,
        passAtmosphere=True,
        reuse_wavelength_grid=transfer_grid == "fixed",
    )

    evaluation_sec = []
    evaluation_flux = []
    for offset in offsets:
        sme.abund.A["Fe"] = base_fe + offset
        kwargs = {
            "passLineList": True,
            "passAtmosphere": True,
        }
        if supports_prepared and profile == "default":
            kwargs = {
                "passLineList": False,
                "passAtmosphere": False,
                "passAbund": True,
            }
        kwargs["reuse_wavelength_grid"] = transfer_grid == "fixed"
        wall, flux = synthesize(synth, sme, **kwargs)
        evaluation_sec.append(wall)
        evaluation_flux.append(flux)
        print(
            f"{pysme.__version__} {profile} {model[0]} {width:g}A "
            f"Fe={offset:+.2f} {wall:.6f}s",
            flush=True,
        )

    if warm_control is not None and eos_warm:
        warm_control(False)

    values = np.asarray(evaluation_sec, dtype=float)
    return {
        "input_lines": int(len(linelist)),
        "wave_points": int(initial_flux.size),
        "initial_synthesis_sec": initial_sec,
        "evaluation_sec": values.tolist(),
        "median_evaluation_sec": (
            float(np.median(values)) if values.size else None
        ),
        "total_evaluation_sec": float(np.sum(values)),
        "total_workload_sec": float(initial_sec + np.sum(values)),
        "prepared_state_used": bool(supports_prepared and profile == "default"),
        "initial_flux": initial_flux,
        "evaluation_flux": evaluation_flux,
    }


def run(args: argparse.Namespace):
    output = {
        "version": str(pysme.__version__),
        "profile": args.profile,
        "center_A": CENTER_A,
        "step_A": STEP_A,
        "transfer_grid": args.transfer_grid,
        "fe_offsets_dex": list(args.offsets),
        "models": {},
    }
    spectra = {}
    for width in args.widths:
        path = window_path(width)
        if not path.exists():
            raise FileNotFoundError(path)
        started = time.perf_counter()
        linelist = ValdFile(str(path))
        parse_sec = time.perf_counter() - started
        print(
            f"loaded {width:g}A: {len(linelist):,} lines in {parse_sec:.2f}s",
            flush=True,
        )
        selected_models = [model for model in MODELS if model[0] in args.models]
        for model in selected_models:
            record = run_case(
                model,
                linelist,
                width,
                args.profile,
                tuple(args.offsets),
                args.eos_warm,
                args.transfer_grid,
            )
            initial_flux = record.pop("initial_flux")
            evaluation_flux = record.pop("evaluation_flux")
            record["parse_sec_excluded"] = parse_sec
            key = f"{model[0]}_{width:g}"
            output["models"][key] = record
            spectra[f"{key}_initial"] = initial_flux
            for index, flux in enumerate(evaluation_flux):
                spectra[f"{key}_eval_{index}"] = flux
            args.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
        del linelist

    args.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
    np.savez_compressed(args.spectra, **spectra)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--profile", choices=("default", "reference"), default="default")
    parser.add_argument(
        "--transfer-grid",
        choices=("fixed", "wave-only"),
        default="fixed",
    )
    parser.add_argument("--widths", type=float, nargs="+", default=[10.0, 200.0])
    parser.add_argument(
        "--offsets", type=float, nargs="*", default=list(FE_OFFSETS_DEX)
    )
    parser.add_argument(
        "--models", nargs="+", default=[model[0] for model in MODELS]
    )
    parser.add_argument(
        "--eos-warm", action=argparse.BooleanOptionalAction, default=True
    )
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--spectra", type=Path, required=True)
    return parser.parse_args()


if __name__ == "__main__":
    run(parse_args())
