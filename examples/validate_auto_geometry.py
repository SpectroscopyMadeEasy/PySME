"""Exercise automatic MARCS geometry through synthesis and synthetic fits.

Run from a source checkout after installing it:
    python examples/validate_auto_geometry.py --out /path/to/new/results

The output compares the actual automatic synthesis pipeline across logg=3.5
and fits noiseless synthetic spectra from initial gravities on both sides.
It is a numerical consistency diagnostic, not an observational validation.
"""
import argparse
from copy import deepcopy
import hashlib
import json
from pathlib import Path
import subprocess

import numpy as np
import pandas as pd

from pysme import util
from pysme.abund import Abund
from pysme.linelist.linelist import LineList
from pysme.linelist.vald import ValdFile
from pysme.sme import SME_Structure
from pysme.solve import SME_Solver
from pysme.synthesize import Synthesizer


ROOT = Path(__file__).resolve().parents[1]
CASES = [(4000., .5), (5500., 0.), (7000., -2.5)]
TRUTHS = [3.45, 3.499, 3.5, 3.501, 3.55]
WINDOWS = [(6436., 6444.), (6550., 6575.)]


class TracingSynthesizer(Synthesizer):
    def __init__(self):
        super().__init__()
        self.trace = []
        self.limit = None

    def get_atmosphere(self, sme):
        if self.limit is not None and len(self.trace) >= self.limit:
            raise RuntimeError("Diagnostic fit exceeded 100 atmosphere evaluations")
        result = super().get_atmosphere(sme)
        assert result.atmo.geom_request is None
        expected = "SPH" if result.logg < 3.5 else "PP"
        assert result.atmo.geom == expected
        self.trace.append({"logg": float(result.logg), "geom": result.atmo.geom})
        return result


def make_sme(teff, monh, source, lines, nmu):
    sme = SME_Structure()
    sme.teff, sme.monh = teff, monh
    sme.abund = Abund(monh=monh, pattern="grevesse2007")
    sme.vmic, sme.vmac, sme.vsini = 2., 0., 0.
    sme.atmo.source, sme.atmo.geom = source, None
    sme.atmo.depth, sme.atmo.interp = "RHOX", "RHOX"
    sme.linelist = deepcopy(lines)
    sme.wave = [np.round(np.arange(a, b + .005, .01), 6) for a, b in WINDOWS]
    sme.mu = np.sqrt((2*np.arange(nmu)+1)/(2*nmu))
    sme.line_select_method = "internal"
    sme.accrt, sme.accwi = 1e-5, 1e-5
    sme.vrad_flag, sme.cscale_flag = "none", "none"
    sme.iptype, sme.ipres = "gauss", 50000.
    sme.accxt, sme.accft, sme.accgt = 1e-5, 1e-8, 1e-8
    return sme


def profiles(sme):
    return [np.asarray(f).copy() for f in sme.synth]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True, help="New output directory")
    parser.add_argument("--source", default="marcs2014.sav")
    parser.add_argument("--mu", type=int, default=7)
    parser.add_argument("--case", type=int, choices=range(len(CASES)))
    parser.add_argument("--skip-fits", action="store_true")
    args = parser.parse_args()
    if not 1 <= args.mu <= 77:
        parser.error("--mu must be between 1 and the native limit of 77")
    args.out.mkdir(parents=True, exist_ok=False)
    util.show_progress_bars = False
    paths = [ROOT/"test/testcase1.lin", ROOT/"test/testcase3.lin"]
    frames = [ValdFile(str(p))._lines for p in paths]
    lines = LineList(pd.concat(frames, ignore_index=True).sort_values("wlcent").reset_index(drop=True),
                     lineformat="short", medium="air")
    synth = TracingSynthesizer()
    cases = CASES if args.case is None else [CASES[args.case]]
    gravities = np.unique(np.r_[np.linspace(3.45, 3.55, 51), TRUTHS,
                               3.5-1e-6, 3.5, 3.5+1e-6])
    results, fits = [], []
    all_traces = {}
    for teff, monh in cases:
        label = f"T{teff:.0f}_M{monh:+.2f}"
        sme = make_sme(teff, monh, args.source, lines, args.mu)
        spectra = {}
        for gravity in gravities:
            sme.logg = float(gravity)
            synth.synthesize_spectrum(sme)
            spectra[float(gravity)] = profiles(sme)
        # Repeated switches exercise state reuse rather than independent models.
        for gravity in [3.51, 3.49, 3.51, 3.49]:
            sme.logg = gravity
            synth.synthesize_spectrum(sme)
            reference_key = min(spectra, key=lambda key: abs(key-gravity))
            for actual, expected in zip(profiles(sme), spectra[reference_key]):
                np.testing.assert_allclose(actual, expected, atol=1e-9, rtol=1e-9)
        for gravity in gravities:
            row = {"case":label, "teff":teff, "monh":monh, "logg":float(gravity),
                   "geom":"SPH" if gravity < 3.5 else "PP"}
            for truth in [3.499, 3.501]:
                residual = np.concatenate([a-b for a,b in zip(spectra[float(gravity)], spectra[truth])])
                row[f"chi2_truth_{truth}"] = float(np.sum((residual*300)**2))
            results.append(row)
        data = {"logg":gravities}
        for i in range(2):
            data[f"wave{i}"] = np.asarray(sme.wave[i])
            data[f"flux{i}"] = np.array([spectra[float(g)][i] for g in gravities])
        np.savez_compressed(args.out/f"{label}.npz", **data)
        print(f"{label}: {len(gravities)} spectra; repeated switches match", flush=True)
        if not args.skip_fits:
            for truth in TRUTHS:
                for initial in [3.3, 3.7]:
                    trial = make_sme(teff, monh, args.source, lines, args.mu)
                    trial.spec = [f.copy() for f in spectra[truth]]
                    trial.uncs = [np.full_like(f, 1/300) for f in spectra[truth]]
                    trial.logg = initial
                    solver = SME_Solver()
                    tracer = TracingSynthesizer()
                    tracer.limit = 100
                    solver.synthesizer = tracer
                    fit = {"case":label,"truth":truth,"initial":initial}
                    try:
                        solved = solver.solve(trial, ["logg"], bounds=np.array([[3.2],[3.8]]))
                        assert solved.atmo.geom == ("SPH" if solved.logg < 3.5 else "PP")
                        assert solved.atmo.geom_request is None
                        fit.update(logg=float(solved.logg), delta_logg=float(solved.logg-truth),
                                   chisq=float(solved.fitresults.chisq), geom=solved.atmo.geom,
                                   status="completed")
                    except Exception as error:
                        fit.update(status="failed",error=repr(error))
                    fit["evaluations"] = len(tracer.trace)
                    fit["visited_geometries"] = sorted({r["geom"] for r in tracer.trace})
                    fits.append(fit)
                    all_traces[f"{label}_truth{truth}_initial{initial}"] = tracer.trace
                    print(json.dumps(fit), flush=True)
                    (args.out/"fit_traces.json").write_text(json.dumps(all_traces,indent=2)+"\n")
                    pd.DataFrame(fits).to_csv(args.out/"fits.csv",index=False)
    pd.DataFrame(results).to_csv(args.out/"curve.csv",index=False)
    try:
        commit = subprocess.check_output(["git","-C",str(ROOT),"rev-parse","HEAD"],text=True).strip()
    except (OSError, subprocess.CalledProcessError):
        commit = "unknown"
    code_paths = [ROOT / p for p in (
        "src/pysme/atmosphere/atmosphere.py", "src/pysme/atmosphere/interpolation.py",
        "src/pysme/synthesize.py", "src/pysme/solve.py", "examples/validate_auto_geometry.py")]
    source_path = Path(args.source)
    metadata = {"commit":commit,"source":args.source,"nmu":args.mu,"R":50000,
                "source_sha256":hashlib.sha256(source_path.read_bytes()).hexdigest() if source_path.is_file() else None,
                "code_sha256":{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in code_paths},
                "cases":cases,"truths":TRUTHS,"initial_logg":[3.3,3.7],"LTE":True,
                "accrt":1e-5,"accwi":1e-5,"jacobian_step_logg":0.01,
                "noise":"none; uncs=1/300 weights only; logg is the only free parameter",
                "linelist_sha256":{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
                "native_version":synth.get_dll().SMELibraryVersion(),
                "failed_fits":sum(f["status"]!="completed" for f in fits),
                "fits_with_abs_delta_logg_above_0.003":sum(abs(f.get("delta_logg",np.inf))>.003 for f in fits)}
    (args.out/"metadata.json").write_text(json.dumps(metadata,indent=2)+"\n")
    print(json.dumps(metadata,indent=2),flush=True)


if __name__ == "__main__":
    main()
