"""Benchmark a coarse wavelength cache for SMElib continuum opacity.

This is analysis-only code.  The CDR/spectrum modes expect the transient native
BenchmarkContinuum* and BenchmarkGetLineALMAX instrumentation described in the
accompanying report.  That instrumentation was removed from the production
source after the saved benchmark artifacts were generated; continuum mode can
still run against the clean library.
"""

from __future__ import annotations

import argparse
import ctypes
import json
import os
import time
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from pysme import util
from pysme.abund import Abund
from pysme.linelist.vald import ValdFile
from pysme.sme import SME_Structure
from pysme.sme_synth import SME_DLL
from pysme.smelib.libtools import get_full_libfile
from pysme.synthesize import Synthesizer


ROOT = Path(__file__).resolve().parents[1]
OUTDIR = ROOT / "analysis"
FULL_LINELIST = Path("/tmp/pysme_cdr_5045_5355.lin")
SMALL_LINELIST = ROOT / "test" / "testcase1.lin"
CONTINUUM_JSON = OUTDIR / "continuum_opacity_grid_continuum.json"
CONTINUUM_NPZ = OUTDIR / "continuum_opacity_grid_reference.npz"
CONTINUUM_PNG = OUTDIR / "continuum_opacity_grid_worst.png"
CDR_JSON = OUTDIR / "continuum_opacity_grid_cdr.json"
CDR_NPZ = OUTDIR / "continuum_opacity_grid_cdr_reference.npz"
SPECTRUM_JSON = OUTDIR / "continuum_opacity_grid_spectrum_windows.json"
SPECTRUM_NPZ = OUTDIR / "continuum_opacity_grid_spectrum_windows.npz"

MODELS = [
    ("solar_dwarf", 5777.0, 4.44, 0.0, 1.0),
    ("metal_poor_dwarf", 6000.0, 4.0, -2.0, 1.2),
    ("cool_metal_rich_dwarf", 4250.0, 4.5, 0.3, 1.0),
    ("giant", 4500.0, 2.0, 0.0, 1.5),
]
STEPS = (2.0, 1.0, 0.5, 0.25)
INTERPS = ("linear", "log")

# The first two are true H I bound-free thresholds computed from the constants
# used by COULX.  The remaining pairs are sharp/high-curvature PEACH table
# knots in the Mg I and Si I implementations.
DISCONTINUITIES = np.array(
    [
        3647.055245510257,
        8205.874302398079,
        3756.6077790090994,
        6549.676211125072,
        7234.20444455255,
        7291.823802752545,
        5379.96003045514,
        5625.05598459199,
        6260.486002117759,
        6349.3269354463155,
        6491.10396629149,
    ],
    dtype=float,
)
FORCED_NODES = np.sort(
    np.concatenate(
        [DISCONTINUITIES, DISCONTINUITIES + 1e-4,
         np.array([3756., 3757., 6549., 6550., 7234., 7235., 7291., 7292.,
                   5379., 5380., 5624., 5625., 6260., 6261., 6349., 6350.,
                   6491., 6492.])]
    )
)
LOCAL_EDGE_OFFSETS = np.array(
    [-1.0, -0.5, -0.25, -0.1, -0.05, -0.02, 0.0,
     1e-4, 0.02, 0.05, 0.1, 0.25, 0.5, 1.0],
    dtype=float,
)

WAVE = np.linspace(5195.0, 5205.0, 201)
CDR_WAVE = np.arange(5000.0, 5010.0, 1.0)
SPECTRUM_WINDOWS = {
    "blue_3870_3880": (Path("/tmp/pysme_padded_vald/window_3870_3880_full.lin"), 3870.0, 3880.0),
    "blue_4495_4505": (Path("/tmp/pysme_padded_vald/window_4495_4505_full.lin"), 4495.0, 4505.0),
    "optical_5145_5155": (Path("/tmp/pysme_padded_vald/window_5145_5155_full.lin"), 5145.0, 5155.0),
    "nad_5885_5895": (Path("/tmp/pysme_nad_vald/window_5885_5895_full.lin"), 5885.0, 5895.0),
    "halpha_6560_6570": (Path("/tmp/pysme_padded_vald/window_6560_6570_full.lin"), 6560.0, 6570.0),
    "red_7545_7555": (Path("/tmp/pysme_padded_vald/window_7545_7555_full.lin"), 7545.0, 7555.0),
    "caii_8545_8555": (Path("/tmp/pysme_padded_vald/window_8545_8555_full.lin"), 8545.0, 8555.0),
}

util.show_progress_bars = False


def native_api():
    lib = ctypes.CDLL(get_full_libfile())
    lib.BenchmarkContinuumCacheReset.argtypes = []
    lib.BenchmarkContinuumCacheReset.restype = None
    lib.BenchmarkContinuumCacheStats.argtypes = [
        ctypes.POINTER(ctypes.c_longlong),
        ctypes.POINTER(ctypes.c_longlong),
        ctypes.POINTER(ctypes.c_longlong),
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double),
    ]
    lib.BenchmarkContinuumCacheStats.restype = None
    lib.BenchmarkGetLineALMAX.argtypes = [
        ctypes.c_int,
        np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags="C_CONTIGUOUS"),
    ]
    lib.BenchmarkGetLineALMAX.restype = ctypes.c_int
    return lib


def reset_cache(lib):
    lib.BenchmarkContinuumCacheReset()


def cache_stats(lib):
    values = [ctypes.c_longlong() for _ in range(3)]
    times = [ctypes.c_double() for _ in range(2)]
    lib.BenchmarkContinuumCacheStats(
        *(ctypes.byref(v) for v in values), *(ctypes.byref(v) for v in times)
    )
    return {
        "queries": values[0].value,
        "exact_calls": values[1].value,
        "resident_nodes": values[2].value,
        "grid_construction_sec": times[0].value,
        "interpolation_sec": times[1].value,
    }


def set_cache(step=0.0, interp="linear", edge_aware=False):
    os.environ["SME_CONT_CACHE_STEP"] = str(float(step))
    os.environ["SME_CONT_CACHE_INTERP"] = interp
    os.environ["SME_CONT_CACHE_EDGE_AWARE"] = "1" if edge_aware else "0"


def make_sme(name, teff, logg, monh, vmic, linelist):
    sme = SME_Structure()
    sme.teff = teff
    sme.logg = logg
    sme.monh = monh
    sme.vmic = vmic
    sme.vmac = 0.0
    sme.vsini = 0.0
    sme.ipres = 1.0e7
    sme.abund = Abund(monh=monh, pattern="asplund2009")
    sme.linelist = ValdFile(str(linelist))
    sme.atmo.source = "marcs2012.sav"
    sme.atmo.method = "grid"
    sme.wave = [WAVE]
    sme.wint = [WAVE]
    sme.vrad_flag = "none"
    sme.cscale_flag = "none"
    sme.accrt = 1e-4
    sme.line_select_method = "cdr"
    sme.line_select_policy = "strict"
    sme.line_select_parallel = False
    sme.line_select_recompute = "never"
    sme.line_select_cdr_strength_thres = 0.001
    sme.line_select_cdr_bin_width = 0.2
    sme.strong_depth_thres = 0.001
    sme.strong_bin_width = 0.2
    return Synthesizer().get_atmosphere(sme)


def prepare_dll(sme, wave_range):
    dll = SME_DLL()
    dll.SetLineInfoMode(0)
    ion_mask = np.asarray(dll.InputLineList(sme.linelist), dtype=bool)
    dll.InputModel(sme.teff, sme.logg, sme.vmic, sme.atmo)
    dll.InputAbund(sme.abund)
    t0 = time.perf_counter()
    dll.Ionization(0)
    ionization_sec = time.perf_counter() - t0
    dll.SetVWscale(sme.gam6)
    dll.SetH2broad(sme.h2broad)
    dll.InputWaveRange(*wave_range)
    dll.Opacity()
    return dll, ion_mask, ionization_sec


def grid_nodes(step, edge_aware, lo=3500.0, hi=9500.0):
    nodes = np.arange(np.floor(lo / step) * step, hi + step * 1.01, step)
    if edge_aware:
        extra = FORCED_NODES[(FORCED_NODES >= lo) & (FORCED_NODES <= hi)]
        nodes = np.unique(np.concatenate([nodes, extra]))
    return nodes


def interpolate(nodes, values, wave, method, edge_aware=False):
    right = np.searchsorted(nodes, wave, side="left")
    right = np.clip(right, 1, len(nodes) - 1)
    left = right - 1
    exact_right = np.isclose(wave, nodes[right], rtol=0.0, atol=1e-12)
    left[exact_right] = right[exact_right]
    if edge_aware:
        for edge in DISCONTINUITIES:
            blue = int(np.searchsorted(nodes, edge, side="left"))
            red = blue + 1
            just_red = (wave > edge) & (wave < nodes[red])
            left[just_red] = red
            right[just_red] = min(red + 1, len(nodes) - 1)
    denom = nodes[right] - nodes[left]
    frac = np.zeros_like(wave)
    nz = denom != 0
    frac[nz] = (wave[nz] - nodes[left[nz]]) / denom[nz]
    a = values[left]
    b = values[right]
    f = frac[:, None]
    if method == "log":
        positive = (a > 0) & (b > 0)
        linear = a + f * (b - a)
        logged = np.exp(np.log(np.maximum(a, 1e-300)) + f * (np.log(np.maximum(b, 1e-300)) - np.log(np.maximum(a, 1e-300))))
        return np.where(positive, logged, linear)
    return a + f * (b - a)


def error_summary(approx, exact, wave):
    absolute = np.abs(approx - exact)
    scale = np.maximum(np.abs(exact), np.nanmedian(np.abs(exact), axis=0) * 1e-12)
    relative = absolute / np.maximum(scale, 1e-300)
    flat = int(np.nanargmax(relative))
    iw, depth = np.unravel_index(flat, relative.shape)
    return {
        "absolute": {
            "median": float(np.nanmedian(absolute)),
            "p90": float(np.nanpercentile(absolute, 90)),
            "p99": float(np.nanpercentile(absolute, 99)),
            "max": float(np.nanmax(absolute)),
        },
        "relative": {
            "median": float(np.nanmedian(relative)),
            "p90": float(np.nanpercentile(relative, 90)),
            "p99": float(np.nanpercentile(relative, 99)),
            "max": float(relative[iw, depth]),
        },
        "worst_wavelength_A": float(wave[iw]),
        "worst_depth_index": int(depth),
    }


def band_error_summary(approx, exact, wave):
    bands = {
        "blue_3500_4500": (3500.0, 4500.0),
        "visible_4500_7000": (4500.0, 7000.0),
        "red_7000_9500": (7000.0, 9500.0),
    }
    return {
        name: error_summary(approx[(wave >= lo) & (wave < hi)],
                            exact[(wave >= lo) & (wave < hi)],
                            wave[(wave >= lo) & (wave < hi)])
        for name, (lo, hi) in bands.items()
    }


def continuum_eval_wave():
    regular = np.arange(3500.125, 9500.0, 0.5)
    dense = [np.arange(edge - 1.0, edge + 1.0001, 0.02) for edge in FORCED_NODES]
    return np.unique(np.concatenate([regular, *dense]))


def run_continuum():
    output = {"models": {}, "steps_A": STEPS, "interpolation": INTERPS}
    reference_npz = {}
    plot_rows = []
    eval_wave = continuum_eval_wave()
    for model in MODELS:
        name, teff, logg, monh, vmic = model
        print(f"continuum: {name}", flush=True)
        sme = make_sme(*model, SMALL_LINELIST)
        dll, _, _ = prepare_dll(sme, (3490.0, 9510.0))
        memo = {}

        def exact_at(waves):
            rows = []
            for wavelength in waves:
                key = float(wavelength)
                if key not in memo:
                    memo[key] = tuple(
                        np.asarray(x, dtype=float)
                        for x in dll.GetContinuumOpacityComponents(key)
                    )
                rows.append(memo[key])
            return tuple(np.stack([row[i] for row in rows]) for i in range(3))

        t0 = time.perf_counter()
        kappa_ref, sigma_ref, chi_ref = exact_at(eval_wave)
        model_out = {
            "atmosphere": {"teff": teff, "logg": logg, "monh": monh, "vmic": vmic},
            "n_eval_wavelength": len(eval_wave),
            "n_depth": int(chi_ref.shape[1]),
            "configs": {},
        }
        output["models"][name] = model_out
        for step in STEPS:
            for edge_aware in (False, True):
                nodes = grid_nodes(step, edge_aware)
                kappa_nodes, sigma_nodes, chi_nodes = exact_at(nodes)
                for interp in INTERPS:
                    key = f"{step:g}A_{interp}_{'edge' if edge_aware else 'naive'}"
                    kappa = interpolate(nodes, kappa_nodes, eval_wave, interp, edge_aware)
                    sigma = interpolate(nodes, sigma_nodes, eval_wave, interp, edge_aware)
                    chi = interpolate(nodes, chi_nodes, eval_wave, interp, edge_aware)
                    model_out["configs"][key] = {
                        "nodes": len(nodes),
                        "kappa": error_summary(kappa, kappa_ref, eval_wave),
                        "sigma": error_summary(sigma, sigma_ref, eval_wave),
                        "chi": error_summary(chi, chi_ref, eval_wave),
                        "chi_bands": band_error_summary(chi, chi_ref, eval_wave),
                    }
                    if name == "solar_dwarf":
                        per_wave = np.max(np.abs((chi - chi_ref) / np.maximum(np.abs(chi_ref), 1e-300)), axis=1)
                        plot_rows.append((key, per_wave))
        refined_nodes = np.unique(
            np.concatenate(
                [grid_nodes(1.0, True),
                 (DISCONTINUITIES[:, None] + LOCAL_EDGE_OFFSETS).ravel()]
            )
        )
        refined_nodes = refined_nodes[(refined_nodes >= 3500.0) & (refined_nodes <= 9500.0)]
        kappa_nodes, sigma_nodes, chi_nodes = exact_at(refined_nodes)
        for interp in INTERPS:
            key = f"1A_{interp}_locally_refined"
            kappa = interpolate(refined_nodes, kappa_nodes, eval_wave, interp, True)
            sigma = interpolate(refined_nodes, sigma_nodes, eval_wave, interp, True)
            chi = interpolate(refined_nodes, chi_nodes, eval_wave, interp, True)
            model_out["configs"][key] = {
                "nodes": len(refined_nodes),
                "kappa": error_summary(kappa, kappa_ref, eval_wave),
                "sigma": error_summary(sigma, sigma_ref, eval_wave),
                "chi": error_summary(chi, chi_ref, eval_wave),
                "chi_bands": band_error_summary(chi, chi_ref, eval_wave),
            }
        model_out["exact_wall_sec"] = time.perf_counter() - t0
        model_out["unique_exact_calls"] = len(memo)
        output["models"][name] = model_out
        if name == "solar_dwarf":
            reference_npz.update(
                wave=eval_wave, kappa=kappa_ref, sigma=sigma_ref, chi=chi_ref
            )
        CONTINUUM_JSON.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")

    np.savez_compressed(CONTINUUM_NPZ, **reference_npz)
    fig, ax = plt.subplots(figsize=(12, 6))
    for key, error in plot_rows:
        if key in {"2A_linear_naive", "1A_linear_naive", "1A_log_edge", "0.5A_log_edge", "0.25A_log_edge"}:
            ax.semilogy(eval_wave, np.maximum(error, 1e-16), label=key)
    for edge in FORCED_NODES:
        ax.axvline(edge, color="0.85", linewidth=0.5)
    ax.set(xlabel="Wavelength [Å]", ylabel="max depth relative error in continuum extinction")
    ax.legend(ncol=2)
    fig.tight_layout()
    fig.savefig(CONTINUUM_PNG, dpi=160)
    print(CONTINUUM_JSON)


def calculate_cdr(sme, accrt=1e-4):
    native = native_api()
    dll, ion_mask, ionization_sec = prepare_dll(
        sme, (float(CDR_WAVE[0] - 2), float(CDR_WAVE[-1] + 2))
    )
    reset_cache(native)
    t0 = time.perf_counter()
    dll.Transf(sme.mu, wave=CDR_WAVE, accrt=accrt)
    transf_sec = time.perf_counter() - t0
    ranges_valid = np.asarray(dll.GetLineRange(), dtype=float)
    almax_valid = np.empty(len(ranges_valid), dtype=np.float64)
    copied = native.BenchmarkGetLineALMAX(len(almax_valid), almax_valid)
    if copied != len(almax_valid):
        raise RuntimeError(f"ALMAX copy mismatch {copied} != {len(almax_valid)}")
    t0 = time.perf_counter()
    central_valid = np.asarray(dll.CentralDepth(sme.mu, accrt), dtype=float)
    central_depth_sec = time.perf_counter() - t0
    stats = cache_stats(native)

    keep = ~ion_mask
    n = len(sme.linelist)
    central = np.full(n, np.nan)
    almax = np.full(n, np.nan)
    range_s = np.full(n, np.nan)
    range_e = np.full(n, np.nan)
    central[keep] = central_valid
    almax[keep] = almax_valid
    range_s[keep] = ranges_valid[:, 0]
    range_e[keep] = ranges_valid[:, 1]
    species = np.char.strip(np.asarray(sme.linelist["species"], dtype=str))
    central[species == "H 1"] = 1.0
    widths = range_e - range_s
    saturated = np.isclose(widths, 2000.0, rtol=1e-4, atol=5.0, equal_nan=False)
    wl = np.asarray(sme.linelist["wlcent"], dtype=float)
    range_s[saturated] = wl[saturated] - 0.3
    range_e[saturated] = wl[saturated] + 0.3
    strong = np.asarray(
        Synthesizer.flag_strong_lines_by_bins(
            wl, central, bin_width=0.2, threshold=0.001, valid_mask=keep
        ),
        dtype=bool,
    )
    timing = {
        "ionization_sec": ionization_sec,
        "transf_sec": transf_sec,
        "central_depth_sec": central_depth_sec,
        "total_sec": ionization_sec + transf_sec + central_depth_sec,
        "strong_lines": int(strong.sum()),
        "valid_lines": int(keep.sum()),
        "cache": stats,
    }
    return central, almax, range_s, range_e, strong, timing


def install_cdr(sme, central, range_s, range_e, strong):
    frame = sme.linelist._lines
    frame["central_depth"] = central
    frame["line_range_s"] = range_s
    frame["line_range_e"] = range_e
    frame["strong"] = strong
    sme.linelist.cdr_paras = np.array([sme.teff, sme.logg, sme.monh, sme.vmic])
    sme.linelist.cdr_paras_h_stark_convolution = "legacy"
    sme.linelist.cdr_paras_thres["strong_depth"] = 0.001
    sme.linelist.cdr_paras_thres["strong_bin_width"] = 0.2
    sme.accrt = 1e-4


def synthesize(sme, arrays, step, interp, edge_aware):
    central, _, range_s, range_e, strong = arrays
    install_cdr(sme, central, range_s, range_e, strong)
    set_cache(step, interp, edge_aware)
    native = native_api()
    reset_cache(native)
    t0 = time.perf_counter()
    Synthesizer().synthesize_spectrum(
        sme, linelist_mode="dynamic", smelib_lineinfo_mode=2
    )
    elapsed = time.perf_counter() - t0
    return np.asarray(sme.synth[0], dtype=float).copy(), {
        "sec": elapsed,
        "selected_lines": int(np.count_nonzero(sme.linelist["use_indices"])),
        "cache": cache_stats(native),
    }


def vector_error(approx, exact, relative=True):
    mask = np.isfinite(approx) & np.isfinite(exact)
    delta = np.abs(approx[mask] - exact[mask])
    if relative:
        delta = delta / np.maximum(np.abs(exact[mask]), 1e-300)
    return {
        "median": float(np.median(delta)),
        "p90": float(np.percentile(delta, 90)),
        "p99": float(np.percentile(delta, 99)),
        "max": float(np.max(delta)),
    }


def spectrum_error(approx, exact, wave=WAVE):
    delta = approx - exact
    i = int(np.argmax(np.abs(delta)))
    return {
        "max_abs": float(np.max(np.abs(delta))),
        "rms": float(np.sqrt(np.mean(delta * delta))),
        "p99_abs": float(np.percentile(np.abs(delta), 99)),
        "mean": float(np.mean(delta)),
        "worst_wavelength_A": float(wave[i]),
    }


def selection_comparison(approx, exact):
    return {
        "exact_selected": int(exact.sum()),
        "approx_selected": int(approx.sum()),
        "false_positive": int(np.count_nonzero(approx & ~exact)),
        "false_negative": int(np.count_nonzero(exact & ~approx)),
        "fraction_changed": float(np.count_nonzero(approx != exact) / len(exact)),
    }


def threshold_refinement(sme, approx_central, exact_central, exact_strong):
    wl = np.asarray(sme.linelist["wlcent"], dtype=float)
    valid = np.isfinite(approx_central) & np.isfinite(exact_central)
    out = {}
    for margin in (0.1, 0.2, 0.3):
        lo, hi = 0.001 * (1 - margin), 0.001 * (1 + margin)
        refine = valid & (approx_central >= lo) & (approx_central <= hi)
        metric = approx_central.copy()
        metric[refine] = exact_central[refine]
        selected = np.asarray(
            Synthesizer.flag_strong_lines_by_bins(
                wl, metric, bin_width=0.2, threshold=0.001, valid_mask=valid
            ),
            dtype=bool,
        )
        out[f"{int(margin * 100)}pct"] = {
            "exact_recomputations": int(refine.sum()),
            **selection_comparison(selected, exact_strong),
        }
    return out


def run_cdr():
    if not FULL_LINELIST.exists():
        raise FileNotFoundError(FULL_LINELIST)
    output = {"linelist": str(FULL_LINELIST), "models": {}}
    reference_npz = {"wave": WAVE}
    for model in MODELS:
        name = model[0]
        print(f"cdr exact: {name}", flush=True)
        set_cache(0.0)
        sme_ref = make_sme(*model, FULL_LINELIST)
        exact = calculate_cdr(sme_ref)
        ref_spec, ref_synth = synthesize(sme_ref, exact[:5], 0.0, "linear", False)
        model_out = {
            "atmosphere": dict(zip(("name", "teff", "logg", "monh", "vmic"), model)),
            "exact": exact[5],
            "exact_synthesis": ref_synth,
            "configs": {},
        }
        output["models"][name] = model_out
        CDR_JSON.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
        reference_npz[f"{name}_central"] = exact[0]
        reference_npz[f"{name}_almax"] = exact[1]
        reference_npz[f"{name}_range_s"] = exact[2]
        reference_npz[f"{name}_range_e"] = exact[3]
        reference_npz[f"{name}_strong"] = exact[4]
        reference_npz[f"{name}_spectrum"] = ref_spec

        configs = [(step, interp, True) for step in STEPS for interp in INTERPS]
        if name == "solar_dwarf":
            configs += [(2.0, interp, False) for interp in INTERPS]
            configs += [(1.0, interp, False) for interp in INTERPS]
        for step, interp, edge_aware in configs:
            key = f"{step:g}A_{interp}_{'edge' if edge_aware else 'naive'}"
            print(f"cdr {name}: {key}", flush=True)
            set_cache(step, interp, edge_aware)
            sme = make_sme(*model, FULL_LINELIST)
            approx = calculate_cdr(sme)
            total_spec, synth_stats = synthesize(
                sme, approx[:5], step, interp, edge_aware
            )
            exact_synth_spec, exact_synth_stats = synthesize(
                sme, approx[:5], 0.0, "linear", False
            )
            changed = approx[4] != exact[4]
            species = np.asarray(sme.linelist["species"], dtype=str)
            wl = np.asarray(sme.linelist["wlcent"], dtype=float)
            worst = np.argsort(
                np.nan_to_num(
                    np.abs(approx[0] - exact[0]) / np.maximum(np.abs(exact[0]), 1e-300),
                    nan=-1,
                )
            )[-10:][::-1]
            model_out["configs"][key] = {
                "timing": approx[5],
                "central_depth_relative": vector_error(approx[0], exact[0]),
                "almax_relative": vector_error(approx[1], exact[1]),
                "range_left_abs_A": vector_error(approx[2], exact[2], relative=False),
                "range_right_abs_A": vector_error(approx[3], exact[3], relative=False),
                "range_changed": int(
                    np.count_nonzero(
                        np.isfinite(approx[2])
                        & np.isfinite(exact[2])
                        & ((approx[2] != exact[2]) | (approx[3] != exact[3]))
                    )
                ),
                "selection": selection_comparison(approx[4], exact[4]),
                "selection_changed_species": {
                    str(item): int(count)
                    for item, count in zip(*np.unique(species[changed], return_counts=True))
                },
                "spectrum_total": spectrum_error(total_spec, ref_spec),
                "spectrum_selection_only": spectrum_error(exact_synth_spec, ref_spec),
                "synthesis": synth_stats,
                "exact_synthesis_with_approx_selection": exact_synth_stats,
                "threshold_refinement": threshold_refinement(
                    sme, approx[0], exact[0], exact[4]
                ),
                "worst_lines": [
                    {
                        "index": int(i),
                        "wavelength_A": float(wl[i]),
                        "species": str(species[i]).strip(),
                        "exact_central_depth": float(exact[0][i]),
                        "approx_central_depth": float(approx[0][i]),
                    }
                    for i in worst
                    if np.isfinite(exact[0][i]) and np.isfinite(approx[0][i])
                ],
            }
            CDR_JSON.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
        CDR_JSON.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
        np.savez_compressed(CDR_NPZ, **reference_npz)
    print(CDR_JSON)


def run_spectrum_windows(nad_only=False):
    if nad_only:
        output = json.loads(SPECTRUM_JSON.read_text())
        with np.load(SPECTRUM_NPZ) as saved:
            spectra = {key: saved[key] for key in saved.files}
        models = MODELS[:1]
    else:
        output = {"models": {}, "sampling_A": 0.05}
        spectra = {}
        models = MODELS
    configs = [
        (0.0, "linear", False, "exact"),
        (1.0, "linear", True, "1A_linear_edge"),
        (1.0, "log", True, "1A_log_edge"),
        (0.5, "linear", True, "0.5A_linear_edge"),
    ]
    for model in models:
        name = model[0]
        if nad_only:
            selected_windows = {"nad_5885_5895": SPECTRUM_WINDOWS["nad_5885_5895"]}
        elif name == "solar_dwarf":
            selected_windows = {
                key: value for key, value in SPECTRUM_WINDOWS.items()
                if key != "nad_5885_5895"
            }
        else:
            selected_windows = {
                key: SPECTRUM_WINDOWS[key]
                for key in ("blue_3870_3880", "halpha_6560_6570", "caii_8545_8555")
            }
        model_out = output["models"].setdefault(name, {"windows": {}})
        for window_name, (linelist, lo, hi) in selected_windows.items():
            if not linelist.exists():
                raise FileNotFoundError(linelist)
            print(f"spectrum {name}: {window_name}", flush=True)
            wave = np.linspace(lo, hi, 201)
            window_out = {"linelist": str(linelist), "configs": {}}
            model_out["windows"][window_name] = window_out
            reference = None
            for step, interp, edge_aware, key in configs:
                set_cache(step, interp, edge_aware)
                sme = make_sme(*model, linelist)
                sme.wave = [wave]
                sme.wint = [wave]
                sme.line_select_method = "internal"
                sme.line_select_recompute = "if_stale"
                native = native_api()
                reset_cache(native)
                t0 = time.perf_counter()
                Synthesizer().synthesize_spectrum(
                    sme, linelist_mode="all", smelib_lineinfo_mode=0
                )
                elapsed = time.perf_counter() - t0
                flux = np.asarray(sme.synth[0], dtype=float).copy()
                entry = {"sec": elapsed, "cache": cache_stats(native)}
                spectra[f"{name}_{window_name}_{key}"] = flux
                if reference is None:
                    reference = flux
                    spectra[f"{name}_{window_name}_wave"] = wave
                else:
                    entry["error"] = spectrum_error(flux, reference, wave)
                window_out["configs"][key] = entry
                SPECTRUM_JSON.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
    np.savez_compressed(SPECTRUM_NPZ, **spectra)
    print(SPECTRUM_JSON)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=("continuum", "cdr", "spectrum", "nad", "all"))
    args = parser.parse_args()
    if args.mode in {"continuum", "all"}:
        run_continuum()
    if args.mode in {"cdr", "all"}:
        run_cdr()
    if args.mode in {"spectrum", "all"}:
        run_spectrum_windows()
    if args.mode == "nad":
        run_spectrum_windows(nad_only=True)


if __name__ == "__main__":
    main()
