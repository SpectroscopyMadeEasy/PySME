"""Release-level end-to-end and controlled 2x2 performance benchmark.

The benchmark deliberately separates VALD parsing from the measured
``Synthesizer.synthesize_spectrum`` call.  Every measured run starts with a
fresh SME structure and Synthesizer, while the operating-system file cache is
left warm.  Raw spectra and line-selection state are written as one NPZ per
run so interrupted multi-hour legacy measurements can be resumed.

The stored controlled sequential profiles were measured with a temporary,
non-public audit patch.  It kept a validated precomputed ALMAX mask and Wlim
array immutable while retaining legacy depth-first wavelength scheduling and
full-scan OPMTRX evaluation.  That patch was removed from the release tree
after measurement; historical and production profiles remain reproducible.
"""

from __future__ import annotations

import argparse
import copy
import ctypes
import hashlib
import json
import os
import platform
import re
import subprocess
import tempfile
import threading
import time
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np
import psutil

import pysme
from pysme.abund import Abund
from pysme.atmosphere.marcsfile import MarcsAtmosphere
from pysme.linelist.vald import ValdFile
from pysme.sme import SME_Structure
from pysme.synthesize import Synthesizer


ROOT = Path(__file__).resolve().parents[1]
WINDOW_DIR = Path("/tmp/pysme_line_centric_windows")
SPHERICAL_CACHE = Path("/tmp/pysme_spherical_batched")
MODEL_ROOT = Path("/Users/mingjie/code/TSFitPy/input_files/model_atmospheres/1D")
OUTPUT = ROOT / "analysis" / "v120_release_benchmark.json"
ARRAY_DIR = ROOT / "analysis" / "v120_release_benchmark_arrays"
ACCRT = 1e-4
ACCWI = 3e-3
OUTPUT_STEP_A = 0.05
CONTROLLED_AUDIT_PATCH_AVAILABLE = False


@dataclass(frozen=True)
class Model:
    name: str
    teff: float
    logg: float
    monh: float
    vmic: float
    geometry: str
    atmosphere: str | None = None


@dataclass(frozen=True)
class Case:
    name: str
    model: str
    lo: float
    hi: float
    linelist: str
    workload: str


MODELS = {
    "solar_dwarf": Model("solar_dwarf", 5777.0, 4.44, 0.0, 1.0, "PP"),
    "metal_poor_dwarf": Model("metal_poor_dwarf", 6000.0, 4.0, -2.0, 1.2, "PP"),
    "cool_metal_rich_dwarf": Model(
        "cool_metal_rich_dwarf", 4250.0, 4.5, 0.3, 1.0, "PP"
    ),
    "moderate_giant": Model(
        "moderate_giant",
        4500.0,
        2.0,
        0.0,
        2.0,
        "SPH",
        str(
            MODEL_ROOT
            / "s4500_g+2.0_m1.0_t02_st_z+0.00_a+0.00_c+0.00_n+0.00_o+0.00_r+0.00_s+0.00.mod"
        ),
    ),
    "cool_line_rich_giant": Model(
        "cool_line_rich_giant",
        3500.0,
        1.0,
        0.0,
        2.0,
        "SPH",
        str(
            MODEL_ROOT
            / "s3500_g+1.0_m1.0_t02_st_z+0.00_a+0.00_c+0.00_n+0.00_o+0.00_r+0.00_s+0.00.mod"
        ),
    ),
    "metal_poor_giant": Model(
        "metal_poor_giant",
        4500.0,
        2.0,
        -2.0,
        2.0,
        "SPH",
        str(
            MODEL_ROOT
            / "s4500_g+2.0_m1.0_t02_st_z-2.00_a+0.40_c+0.00_n+0.00_o+0.40_r+0.00_s+0.00.mod"
        ),
    ),
}


LINE_10 = str(WINDOW_DIR / "window_0010A_pad_150A.lin")
LINE_200 = str(WINDOW_DIR / "window_0200A_pad_150A.lin")
LINE_800 = str(WINDOW_DIR / "window_0800A_pad_150A.lin")
LINE_HALPHA = str(SPHERICAL_CACHE / "halpha_6250_6875.lin")


CASES = {
    # Six atmosphere regimes on the canonical abundance-like window.
    **{
        f"{name}_10A": Case(
            f"{name}_10A", name, 5195.0, 5205.0, LINE_10, "10A"
        )
        for name in MODELS
    },
    # Special/strong profiles in both geometries.
    "solar_dwarf_mgb": Case(
        "solar_dwarf_mgb", "solar_dwarf", 5165.0, 5190.0, LINE_10, "Mg_b"
    ),
    "moderate_giant_mgb": Case(
        "moderate_giant_mgb", "moderate_giant", 5165.0, 5190.0, LINE_10, "Mg_b"
    ),
    "solar_dwarf_halpha": Case(
        "solar_dwarf_halpha",
        "solar_dwarf",
        6550.0,
        6575.0,
        LINE_HALPHA,
        "Halpha",
    ),
    "moderate_giant_halpha": Case(
        "moderate_giant_halpha",
        "moderate_giant",
        6550.0,
        6575.0,
        LINE_HALPHA,
        "Halpha",
    ),
    # Workhorse and wide-band scaling checks.
    "solar_dwarf_200A": Case(
        "solar_dwarf_200A", "solar_dwarf", 5100.0, 5300.0, LINE_200, "200A"
    ),
    "cool_metal_rich_dwarf_200A": Case(
        "cool_metal_rich_dwarf_200A",
        "cool_metal_rich_dwarf",
        5100.0,
        5300.0,
        LINE_200,
        "200A",
    ),
    "moderate_giant_200A": Case(
        "moderate_giant_200A",
        "moderate_giant",
        5100.0,
        5300.0,
        LINE_200,
        "200A",
    ),
    "solar_dwarf_800A": Case(
        "solar_dwarf_800A", "solar_dwarf", 4800.0, 5600.0, LINE_800, "800A"
    ),
    "cool_metal_rich_dwarf_800A": Case(
        "cool_metal_rich_dwarf_800A",
        "cool_metal_rich_dwarf",
        4800.0,
        5600.0,
        LINE_800,
        "800A",
    ),
    "moderate_giant_800A": Case(
        "moderate_giant_800A",
        "moderate_giant",
        4800.0,
        5600.0,
        LINE_800,
        "800A",
    ),
    "metal_poor_giant_800A": Case(
        "metal_poor_giant_800A",
        "metal_poor_giant",
        4800.0,
        5600.0,
        LINE_800,
        "800A",
    ),
}


PROFILES = {
    "H": {
        "continuum": "exact",
        "transfer": "legacy",
        "selection": "internal",
        "controlled": False,
    },
    "C0B0": {
        "continuum": "exact",
        "transfer": "legacy",
        "selection": "almax",
        "controlled": True,
    },
    "C1B0": {
        "continuum": "adaptive",
        "transfer": "legacy",
        "selection": "almax",
        "controlled": True,
    },
    "C0B1": {
        "continuum": "exact",
        "transfer": "batched",
        "selection": "almax",
        "controlled": False,
    },
    "C1B1": {
        "continuum": "adaptive",
        "transfer": "batched",
        "selection": "almax",
        "controlled": False,
    },
}


def sha256(path: str | Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as stream:
        for block in iter(lambda: stream.read(8 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def git_output(*args: str, cwd: Path = ROOT) -> str:
    return subprocess.check_output(["git", *args], cwd=cwd, text=True).strip()


def metadata() -> dict:
    cpu = subprocess.run(
        ["sysctl", "-n", "machdep.cpu.brand_string"],
        text=True,
        capture_output=True,
        check=False,
    ).stdout.strip()
    return {
        "timestamp": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        "machine": {
            "cpu": cpu,
            "platform": platform.platform(),
            "logical_cpus": os.cpu_count(),
        },
        "python": platform.python_version(),
        "pysme_reported_version": str(pysme.__version__),
        "numpy": np.__version__,
        "parent_commit": git_output("rev-parse", "HEAD"),
        "parent_diff_sha256": hashlib.sha256(
            subprocess.check_output(["git", "diff", "--binary"], cwd=ROOT)
        ).hexdigest(),
        "smelib_commit": git_output("rev-parse", "HEAD", cwd=ROOT / "smelib"),
        "smelib_diff_sha256": hashlib.sha256(
            subprocess.check_output(
                ["git", "diff", "--binary"], cwd=ROOT / "smelib"
            )
        ).hexdigest(),
        "compiler": {
            "cxx": subprocess.check_output(
                ["g++-16", "--version"], text=True
            ).splitlines()[0],
            "fortran": subprocess.check_output(
                ["gfortran", "--version"], text=True
            ).splitlines()[0],
            "build_flags": "Autotools Release library; g++-16 -O3 plus configured -g -O2",
        },
        "thread_policy": {
            "processes": 1,
            "line_select_parallel": False,
            "OMP_NUM_THREADS": os.environ.get("OMP_NUM_THREADS", "unset"),
        },
        "cache_policy": (
            "VALD parsing excluded; OS file cache warm; fresh SME structure and "
            "Synthesizer for every measured run"
        ),
        "accrt": ACCRT,
        "accwi": ACCWI,
        "continuum_grid": {
            "base_step_A": 1.0,
            "rtol": 1e-3,
            "min_step_A": 1e-3,
        },
    }


def parse_native_timing(text: str) -> dict:
    starts = [m.start() for m in re.finditer(r"SMElib timing \[Transf#", text)]
    block = text[starts[-1] :] if starts else ""
    result = {"raw": block.strip()}
    patterns = {
        "total_sec": r"^  total\s+:\s+([0-9.eE+-]+) s",
        "setup_sec": r"^  setup\s+:\s+([0-9.eE+-]+) s",
        "lineopac_sec": r"^  LINEOPAC\s+:\s+([0-9.eE+-]+) s, calls=(\d+)",
        "range_scan_sec": r"^  range_scan\s+:\s+([0-9.eE+-]+) s, scans=(\d+), lines=(\d+)",
        "rkints_sec": r"^  RKINTS\s+:\s+([0-9.eE+-]+) s, calls=(\d+)",
        "rkints_sph_sec": r"^  RKINTS_sph\s+:\s+([0-9.eE+-]+) s, calls=(\d+)",
        "opmtrx_sec": r"^  OPMTRX\s+:\s+([0-9.eE+-]+) s, calls=(\d+)",
        "tbintg_sec": r"^  TBINTG\s+:\s+([0-9.eE+-]+) s, calls=(\d+)",
        "tbintg_sph_sec": r"^  TBINTG_sph\s+:\s+([0-9.eE+-]+) s, calls=(\d+)",
    }
    for name, pattern in patterns.items():
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
            adaptive_refined=int(match.group(3)),
            adaptive_accepted=int(match.group(4)),
        )
    return result


class StageTimers:
    def __init__(self):
        self.values: dict[str, float] = {}

    def wrap(self, obj, method: str, key: str | None = None):
        original = getattr(obj, method)
        name = key or method

        def measured(*args, **kwargs):
            started = time.perf_counter()
            try:
                return original(*args, **kwargs)
            finally:
                self.values[name] = self.values.get(name, 0.0) + (
                    time.perf_counter() - started
                )

        setattr(obj, method, measured)


class PeakRSS:
    def __enter__(self):
        self.process = psutil.Process()
        self.start = self.process.memory_info().rss
        self.peak = self.start
        self.stop = threading.Event()

        def sample():
            while not self.stop.wait(0.01):
                self.peak = max(self.peak, self.process.memory_info().rss)

        self.thread = threading.Thread(target=sample, daemon=True)
        self.thread.start()
        return self

    def __exit__(self, *_):
        self.stop.set()
        self.thread.join()
        self.peak = max(self.peak, self.process.memory_info().rss)


def build_sme(case: Case, model: Model, linelist: ValdFile, profile: str):
    cfg = PROFILES[profile]
    wave = np.arange(case.lo, case.hi + 0.5 * OUTPUT_STEP_A, OUTPUT_STEP_A)
    sme = SME_Structure()
    sme.teff = model.teff
    sme.logg = model.logg
    sme.monh = model.monh
    sme.vmic = model.vmic
    sme.vmac = 0.0
    sme.vsini = 0.0
    sme.ipres = 60000.0
    sme.abund = Abund(monh=model.monh, pattern="asplund2009")
    sme.linelist = copy.deepcopy(linelist)
    sme.wave = [wave]
    sme.vrad_flag = "none"
    sme.cscale_flag = "none"
    sme.accrt = ACCRT
    sme.accwi = ACCWI
    sme.line_select_method = cfg["selection"]
    sme.line_select_parallel = False
    sme.line_select_recompute = "always"
    sme.continuum_grid = cfg["continuum"]
    sme.continuum_grid_base_step = 1.0
    sme.continuum_grid_rtol = 1e-3
    sme.continuum_grid_min_step = 1e-3
    sme.transfer_grid_method = cfg["transfer"]
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


def run_once(case: Case, profile: str, base_linelist: ValdFile, index: int):
    cfg = PROFILES[profile]
    model = MODELS[case.model]
    if cfg["controlled"]:
        if not CONTROLLED_AUDIT_PATCH_AVAILABLE:
            raise RuntimeError(
                "C0B0/C1B0 require the temporary controlled-baseline audit "
                "patch; use the stored benchmark results"
            )
        os.environ["SME_RKINTS_CONTROLLED_LINE_STATE"] = "1"
    else:
        os.environ.pop("SME_RKINTS_CONTROLLED_LINE_STATE", None)
    os.environ["SME_TIMING"] = "1"

    sme = build_sme(case, model, base_linelist, profile)
    synth = Synthesizer()
    dll = synth.get_dll()
    timers = StageTimers()
    timers.wrap(synth, "get_atmosphere", "atmosphere_interpolation")
    timers.wrap(synth, "integrate_flux", "flux_integration_broadening")
    for method, key in (
        ("InputLineList", "input_linelist"),
        ("InputModel", "input_model"),
        ("InputAbund", "input_abund"),
        ("Ionization", "eos_ionization"),
        ("Opacity", "opacity_preparation"),
        ("ALMAXRange", "almax_range"),
        ("Transf", "transf"),
    ):
        timers.wrap(dll, method, key)

    with tempfile.TemporaryFile(mode="w+b") as stderr_stream:
        saved_stderr = os.dup(2)
        try:
            os.dup2(stderr_stream.fileno(), 2)
            with PeakRSS() as memory:
                started = time.perf_counter()
                synthesis_result = synth.synthesize_spectrum(
                    sme,
                    updateStructure=False,
                    reuse_wavelength_grid=False,
                    radial_velocity_mode="fast",
                    linelist_mode="all",
                    smelib_lineinfo_mode=0,
                )
                total_sec = time.perf_counter() - started
            ctypes.CDLL(None).fflush(None)
        finally:
            os.dup2(saved_stderr, 2)
            os.close(saved_stderr)
        stderr_stream.seek(0)
        timing_text = stderr_stream.read().decode("utf-8", errors="replace")

    stats = dll.GetContinuumOpacityGridStats()
    wave = np.asarray(synthesis_result[0][0], dtype=float)
    flux = np.asarray(synthesis_result[1][0], dtype=float)
    almax = np.asarray(
        sme.linelist["almax_ratio"]
        if "almax_ratio" in sme.linelist._lines.columns
        else np.full(len(sme.linelist), np.nan),
        dtype=float,
    )
    strong = np.asarray(
        sme.linelist["strong"]
        if "strong" in sme.linelist._lines.columns
        else np.zeros(len(sme.linelist), dtype=bool),
        dtype=bool,
    )
    range_s = np.asarray(
        sme.linelist["line_range_s"]
        if "line_range_s" in sme.linelist._lines.columns
        else np.full(len(sme.linelist), np.nan),
        dtype=float,
    )
    range_e = np.asarray(
        sme.linelist["line_range_e"]
        if "line_range_e" in sme.linelist._lines.columns
        else np.full(len(sme.linelist), np.nan),
        dtype=float,
    )
    array_path = ARRAY_DIR / f"{case.name}__{profile}__{index:02d}.npz"
    np.savez_compressed(
        array_path,
        wave=wave,
        flux=flux,
        almax=almax,
        strong=strong,
        range_s=range_s,
        range_e=range_e,
    )
    continuum_payload_bytes = int(stats["nodes"] * len(sme.atmo.rhox) * 13 * 8)
    native = parse_native_timing(timing_text)
    generation_payload_bytes = int(
        native.get("adaptive_probes", 0) * (1 + 2 * len(sme.mu)) * 8
    )
    return {
        "run": index,
        "total_sec": total_sec,
        "stages_sec": timers.values,
        "native_transf": native,
        "continuum_stats": stats,
        "transfer_points": int(native.get("adaptive_probes", 0)),
        "input_lines": int(len(sme.linelist)),
        "selected_lines": int(np.count_nonzero(strong)),
        "peak_rss_bytes": int(memory.peak),
        "rss_increase_bytes": int(max(0, memory.peak - memory.start)),
        "continuum_payload_bytes": continuum_payload_bytes,
        "generation_payload_bytes": generation_payload_bytes,
        "array_file": str(array_path.relative_to(ROOT)),
    }


def summarize_runs(runs: list[dict]) -> dict:
    times = np.asarray([item["total_sec"] for item in runs], dtype=float)
    return {
        "runs": runs,
        "median_sec": float(np.median(times)),
        "min_sec": float(np.min(times)),
        "max_sec": float(np.max(times)),
        "n": int(times.size),
    }


def load_output() -> dict:
    if OUTPUT.exists():
        return json.loads(OUTPUT.read_text())
    return {"metadata": metadata(), "models": {}, "cases": {}}


def save_output(output: dict):
    OUTPUT.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--cases", nargs="+", choices=sorted(CASES), required=True)
    parser.add_argument("--profiles", nargs="+", choices=PROFILES, required=True)
    parser.add_argument("--repeats", type=int, default=1)
    args = parser.parse_args()
    ARRAY_DIR.mkdir(parents=True, exist_ok=True)
    output = load_output()
    output["metadata"] = metadata()

    for case_name in args.cases:
        case = CASES[case_name]
        model = MODELS[case.model]
        line_path = Path(case.linelist)
        if not line_path.exists():
            raise FileNotFoundError(line_path)
        print(f"parsing {case_name}: {line_path}", flush=True)
        parse_started = time.perf_counter()
        base_linelist = ValdFile(str(line_path))
        parse_sec = time.perf_counter() - parse_started
        output["models"][model.name] = asdict(model)
        record = output["cases"].setdefault(
            case_name,
            {
                "definition": asdict(case),
                "geometry": model.geometry,
                "linelist_sha256": sha256(line_path),
                "parse_sec_excluded": parse_sec,
                "profiles": {},
            },
        )
        for profile in args.profiles:
            existing = record["profiles"].get(profile, {}).get("runs", [])
            start_index = len(existing)
            runs = list(existing)
            for offset in range(args.repeats):
                index = start_index + offset
                print(f"running {case_name} {profile} #{index}", flush=True)
                result = run_once(case, profile, base_linelist, index)
                runs.append(result)
                record["profiles"][profile] = summarize_runs(runs)
                save_output(output)
                print(
                    f"finished {case_name} {profile} #{index}: "
                    f"{result['total_sec']:.6f} s",
                    flush=True,
                )
        del base_linelist
    save_output(output)
    print(OUTPUT)


if __name__ == "__main__":
    main()
