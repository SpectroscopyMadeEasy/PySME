#!/usr/bin/env python3
"""Run identical small PySME synthesis and fit probes on every platform."""

from __future__ import annotations

import argparse
import json
import platform
import sys
from pathlib import Path

import numpy as np

from annotations import emit


FIT_WAVE = np.linspace(6439.5, 6441.5, 81)
# Immutable mock observations: PySME 1.2.0/macOS synthesis at [M/H]=+0.2,
# plus a fixed deterministic residual pattern.  These literals, not a
# platform's synthesis, are passed to solve() everywhere.
FIT_SPEC = np.array(
    [
        0.9653653905594179, 0.9690211877937351, 0.972365110481654, 0.9751663673078939, 0.9776696564478931, 0.979876754519261,
        0.9803528982760713, 0.9819785677133895, 0.9836042534939022, 0.985230154496462, 0.9864179221620829, 0.9875631335581169,
        0.9873483556135343, 0.9884243968816123, 0.9891936916219884, 0.9899629779097733, 0.9907322704384101, 0.9901415691998274,
        0.9909108741997165, 0.9916801854370408, 0.9924494999770482, 0.9930251100270056, 0.993529072074827, 0.9926730372323621,
        0.9931770055015877, 0.9936809768790412, 0.9941849513681685, 0.9946889289672766, 0.9951503760908414, 0.994140258305107,
        0.9944901403551196, 0.9948400237133528, 0.9950671866842914, 0.9951935338767985, 0.9939598797264962, 0.9940862025542204,
        0.9938770689373839, 0.9931496232418754, 0.9908228274757416, 0.9821895988267637, 0.9479232420625999, 0.8729798224059091,
        0.7973205062901975, 0.7867299641309548, 0.8507632599811897, 0.9335548663007645, 0.9764204115576246, 0.9897373126585735,
        0.993465870933509, 0.9951258512701862, 0.9960781971257455, 0.9954423141454596, 0.9959116954709918, 0.9962567174646134,
        0.9941881568708882, 0.9819930330905923, 0.9408199343023526, 0.8969770330669977, 0.9049701679971583, 0.9531049257055044,
        0.9879393200333808, 0.9964774321979273, 0.997759020984212, 0.9966595090032827, 0.9968750713955365, 0.9969578654748905,
        0.996922386676403, 0.9967411772299423, 0.9946786052456784, 0.9936850722294001, 0.9901421759537722, 0.9860816496058863,
        0.9814960782871623, 0.9791253070522354, 0.9830717492041257, 0.9881117848113451, 0.9919043312992221, 0.9953622646438107,
        0.9965329609174557, 0.9976592971250089, 0.9969083238288285,
    ],
    dtype=float,
)


def base_structure(repo: Path, linelist: str, *, teff: float, logg: float, monh: float):
    from pysme.abund import Abund
    from pysme.linelist.vald import ValdFile
    from pysme.sme import SME_Structure

    sme = SME_Structure()
    sme.teff = teff
    sme.logg = logg
    sme.monh = monh
    sme.vmic = 1.0 if logg > 3 else 1.7
    sme.vmac = 0.0
    sme.vsini = 0.0
    sme.abund = Abund(monh=monh, pattern="asplund2009")
    sme.atmo.source = "marcs2012.sav"
    sme.atmo.method = "grid"
    sme.linelist = ValdFile(str(repo / linelist))
    sme.vrad_flag = "none"
    sme.cscale_flag = "none"
    sme.normalize_by_continuum = True
    return sme


def synthesize_case(repo: Path, name: str) -> tuple[np.ndarray, np.ndarray, str]:
    from pysme.synthesize import synthesize_spectrum

    if name.startswith("solar"):
        sme = base_structure(repo, "examples/sun.lin", teff=5772, logg=4.44, monh=0.0)
        if name == "solar_fixed":
            sme.wave = [np.linspace(6439.0, 6442.0, 301)]
            kind = "fixed"
        else:
            sme.wran = [[6439.0, 6442.0]]
            kind = "adaptive"
    elif name == "hydrogen_fixed":
        sme = base_structure(repo, "test/testcase1.lin", teff=5772, logg=4.44, monh=0.0)
        sme.wave = [np.linspace(6561.0, 6564.2, 161)]
        kind = "fixed"
    elif name == "cool_giant_fixed":
        sme = base_structure(repo, "test/testcase1.lin", teff=4286, logg=1.64, monh=-0.52)
        sme.wave = [np.linspace(6561.0, 6564.2, 161)]
        kind = "fixed"
    else:
        raise KeyError(name)
    result = synthesize_spectrum(sme)
    return np.asarray(result.wave[0], dtype=float), np.asarray(result.synth[0], dtype=float), kind


def run_fit(repo: Path) -> dict[str, float]:
    from pysme.sme import MASK_VALUES
    from pysme.solve import solve

    sme = base_structure(repo, "examples/sun.lin", teff=5772, logg=4.44, monh=-0.1)
    sme.wave = [FIT_WAVE.copy()]
    sme.spec = [FIT_SPEC.copy()]
    sme.uncs = [np.full(FIT_WAVE.size, 0.002, dtype=float)]
    sme.mask = [np.full(FIT_WAVE.size, MASK_VALUES.LINE, dtype=np.int16)]
    sme.fitresults.maxiter = 30
    result = solve(sme, ["monh"], segments=[0])
    values = np.asarray(result.fitresults.values, dtype=float)
    uncertainty = np.asarray(result.fitresults.fit_uncertainties, dtype=float)
    return {
        "value": float(values[0]),
        "uncertainty": float(uncertainty[0]),
        "chisq": float(result.fitresults.chisq),
        "iterations": int(result.fitresults.iterations),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--platform", default=platform.system().lower())
    parser.add_argument("--only", nargs="*", choices=("solar_fixed", "solar_adaptive", "hydrogen_fixed", "cool_giant_fixed"))
    parser.add_argument("--skip-fit", action="store_true")
    parser.add_argument("--annotations", action="store_true")
    args = parser.parse_args()
    cases = args.only or ["solar_fixed", "solar_adaptive", "hydrogen_fixed", "cool_giant_fixed"]
    payload: dict[str, object] = {
        "platform": np.array(args.platform),
        "python": np.array(sys.version),
        "case_names": np.asarray(cases),
    }
    for name in cases:
        wave, flux, kind = synthesize_case(args.repo.resolve(), name)
        if not np.all(np.isfinite(flux)):
            raise RuntimeError(f"{name} produced non-finite flux")
        payload[f"{name}__wave"] = wave
        payload[f"{name}__flux"] = flux
        payload[f"{name}__grid_kind"] = np.array(kind)
        line = f"case={name} grid={kind} pixels={wave.size} flux_range=[{flux.min():.9g},{flux.max():.9g}]"
        print(line)
        if args.annotations:
            emit("notice", f"Reference {name}", line)
    fit = None if args.skip_fit else run_fit(args.repo.resolve())
    payload["fit_json"] = np.array(json.dumps(fit, sort_keys=True))
    if fit is not None:
        line = f"value={fit['value']:.9g} uncertainty={fit['uncertainty']:.9g} chisq={fit['chisq']:.9g} iterations={fit['iterations']}"
        print(f"fit {line}")
        if args.annotations:
            emit("notice", "Reference fit", line)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(args.output, **payload)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
