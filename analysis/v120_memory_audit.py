"""Stage-resolved RSS audit for wide, line-rich synthesis workloads."""

from __future__ import annotations

import argparse
import copy
import json
import os
import threading
import time
from contextlib import contextmanager
from pathlib import Path

import numpy as np
import psutil

from pysme.abund import Abund
from pysme.atmosphere.marcsfile import MarcsAtmosphere
from pysme.linelist.vald import ValdFile
from pysme.sme import SME_Structure
from pysme.synthesize import Synthesizer, _adaptive_transfer_capacity

from v120_release_benchmark import CASES, MODELS, OUTPUT_STEP_A


ROOT = Path(__file__).resolve().parents[1]


class RSSProfiler:
    def __init__(self):
        self.process = psutil.Process()
        self.current_stage = "idle"
        self.samples: list[tuple[float, str, int]] = []
        self.stage_peaks: dict[str, int] = {}
        self.stop = threading.Event()

    def rss(self) -> int:
        return int(self.process.memory_info().rss)

    def _sample(self):
        while not self.stop.wait(0.01):
            rss = self.rss()
            stage = self.current_stage
            self.samples.append((time.perf_counter(), stage, rss))
            self.stage_peaks[stage] = max(self.stage_peaks.get(stage, 0), rss)

    def start(self):
        self.thread = threading.Thread(target=self._sample, daemon=True)
        self.thread.start()

    def finish(self):
        self.stop.set()
        self.thread.join()

    @contextmanager
    def stage(self, name: str):
        previous = self.current_stage
        self.current_stage = name
        start_rss = self.rss()
        start = time.perf_counter()
        try:
            yield
        finally:
            end_rss = self.rss()
            self.stage_peaks[name] = max(
                self.stage_peaks.get(name, 0), start_rss, end_rss
            )
            self.current_stage = previous

    def wrap(self, obj, method: str, stage: str | None = None):
        original = getattr(obj, method)
        name = stage or method

        def measured(*args, **kwargs):
            with self.stage(name):
                return original(*args, **kwargs)

        setattr(obj, method, measured)


def build_sme(case_name: str, linelist: ValdFile, copy_linelist: bool):
    case = CASES[case_name]
    model = MODELS[case.model]
    sme = SME_Structure()
    sme.teff = model.teff
    sme.logg = model.logg
    sme.monh = model.monh
    sme.vmic = model.vmic
    sme.vmac = 0.0
    sme.vsini = 0.0
    sme.ipres = 60_000.0
    sme.abund = Abund(monh=model.monh, pattern="asplund2009")
    sme.linelist = copy.deepcopy(linelist) if copy_linelist else linelist
    sme.wave = [
        np.arange(case.lo, case.hi + 0.5 * OUTPUT_STEP_A, OUTPUT_STEP_A)
    ]
    sme.vrad_flag = "none"
    sme.cscale_flag = "none"
    sme.accrt = 1e-4
    sme.accwi = 3e-3
    sme.line_select_method = "almax"
    sme.line_select_parallel = False
    sme.line_select_recompute = "always"
    sme.continuum_grid = "adaptive"
    sme.continuum_grid_base_step = 1.0
    sme.continuum_grid_rtol = 1e-3
    sme.continuum_grid_min_step = 1e-3
    sme.transfer_grid_method = "batched"
    if model.geometry == "SPH":
        atmosphere = MarcsAtmosphere(model.atmosphere)
        sme.atmo = atmosphere
        sme.teff = atmosphere.teff
        sme.logg = atmosphere.logg
        sme.monh = model.monh
        sme.vmic = atmosphere.vturb
        sme.abund = copy.deepcopy(atmosphere.abund)
    else:
        sme.atmo.source = "marcs2012.sav"
        sme.atmo.method = "grid"
        sme.atmo.geom = "PP"
    return sme


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--case",
        choices=sorted(CASES),
        required=True,
    )
    parser.add_argument(
        "--linelist-copy", choices=("none", "deepcopy"), default="none"
    )
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--array-output", type=Path)
    args = parser.parse_args()

    case = CASES[args.case]
    profiler = RSSProfiler()
    baseline_rss = profiler.rss()
    profiler.start()
    with profiler.stage("parse_linelist"):
        linelist = ValdFile(case.linelist)
    after_parse_rss = profiler.rss()
    dataframe_shallow = int(
        linelist._lines.memory_usage(index=True, deep=False).sum()
    )
    dataframe_deep = int(linelist._lines.memory_usage(index=True, deep=True).sum())

    with profiler.stage("build_sme"):
        sme = build_sme(
            args.case, linelist, copy_linelist=args.linelist_copy == "deepcopy"
        )
    after_build_rss = profiler.rss()

    synth = Synthesizer()
    dll = synth.get_dll()
    profiler.wrap(synth, "get_atmosphere", "atmosphere_interpolation")
    profiler.wrap(synth, "integrate_flux", "flux_integration_broadening")
    for method in (
        "InputLineList",
        "InputModel",
        "InputAbund",
        "Ionization",
        "Opacity",
        "ALMAXRange",
        "InputLinePrecomputedInfo",
        "Transf",
    ):
        profiler.wrap(dll, method)

    started = time.perf_counter()
    with profiler.stage("synthesize_spectrum"):
        result = synth.synthesize_spectrum(
            sme,
            updateStructure=False,
            reuse_wavelength_grid=False,
            radial_velocity_mode="fast",
            linelist_mode="all",
            smelib_lineinfo_mode=0,
        )
    elapsed = time.perf_counter() - started
    after_synthesis_rss = profiler.rss()
    profiler.finish()

    stats = dll.GetContinuumOpacityGridStats()
    wave = np.asarray(result[0][0])
    native_nodes = len(synth.wint[0]) if 0 in synth.wint else 0
    transfer_capacity = _adaptive_transfer_capacity(
        sme.linelist, case.lo - 2.0, case.hi + 2.0
    )
    nmu = len(sme.mu)
    transfer_output_bytes = transfer_capacity * (1 + 2 * nmu) * 8
    actual_transfer_bytes = native_nodes * (1 + 2 * nmu) * 8
    nrhox = len(sme.atmo.rhox)
    double_line_state_bytes = len(sme.linelist) * nrhox * 3 * 8
    float_line_state_bytes = len(sme.linelist) * nrhox * 3 * 4

    if args.array_output is not None:
        args.array_output.parent.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(
            args.array_output,
            wave=wave,
            flux=np.asarray(result[1][0]),
            almax=np.asarray(sme.linelist["almax_ratio"]),
            strong=np.asarray(sme.linelist["strong"], dtype=bool),
            range_s=np.asarray(sme.linelist["line_range_s"]),
            range_e=np.asarray(sme.linelist["line_range_e"]),
        )

    record = {
        "case": args.case,
        "linelist_copy": args.linelist_copy,
        "elapsed_sec": elapsed,
        "input_lines": len(sme.linelist),
        "depth_layers": nrhox,
        "output_points": len(wave),
        "adaptive_transfer_nodes": native_nodes,
        "estimated_transfer_capacity": transfer_capacity,
        "baseline_rss_bytes": baseline_rss,
        "after_parse_rss_bytes": after_parse_rss,
        "after_build_rss_bytes": after_build_rss,
        "after_synthesis_rss_bytes": after_synthesis_rss,
        "peak_rss_bytes": max(
            max(item[2] for item in profiler.samples),
            max(profiler.stage_peaks.values()),
        ),
        "stage_peak_rss_bytes": profiler.stage_peaks,
        "linelist_dataframe_shallow_bytes": dataframe_shallow,
        "linelist_dataframe_deep_bytes": dataframe_deep,
        "logical_double_lineop_voigt_bytes": double_line_state_bytes,
        "logical_float_lineop_voigt_bytes": float_line_state_bytes,
        "transfer_output_capacity_bytes": transfer_output_bytes,
        "actual_transfer_node_bytes": actual_transfer_bytes,
        "continuum_stats": stats,
        "continuum_payload_bytes": stats["nodes"] * nrhox * 13 * 8,
    }
    args.output.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
    print(json.dumps(record, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
