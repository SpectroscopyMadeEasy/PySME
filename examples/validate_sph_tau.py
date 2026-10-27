"""Compare RHOX/TAU transfer for one controlled spherical atmosphere.

Run from a source checkout: python examples/validate_sph_tau.py --out NEW_DIR
This is a coordinate-consistency diagnostic, not a MARCS accuracy benchmark.
"""
import argparse
import hashlib
import json
from pathlib import Path
import warnings

import numpy as np
from scipy.integrate import cumulative_trapezoid
from scipy.interpolate import PchipInterpolator

from pysme.abund import Abund
from pysme.atmosphere.krzfile import KrzFile
from pysme.linelist.linelist import LineList
from pysme.sme_synth import SME_DLL
from pysme.smelib.libtools import get_full_libfile


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument(
        "--atmosphere", type=Path,
        default=Path(__file__).resolve().parents[1] / "test/testatmo1.krz",
    )
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    warnings.filterwarnings("ignore", message=".*EOS-computed electron.*")
    base = KrzFile(str(args.atmosphere))
    log_mass = np.log(base.rhox)
    dll = SME_DLL()
    if "SPH:TAU" not in dll.GetModelInputCapabilities().split(","):
        raise RuntimeError("Build SMElib with SPH:TAU support before running this diagnostic")
    lines = LineList()
    lines.add("Fe 1", 5502.9931, 0.9582, -3.047, 7.19, -6.22, 239.249)
    lines.add("Cr 2", 5503.5955, 4.1682, -2.117, 8.37, -6.49, 195.248)
    wave = np.linspace(5502.7, 5503.8, 111)
    mu = np.array([1.0, 0.6, 0.2, 0.05])

    def atmosphere(layers):
        model = KrzFile(str(args.atmosphere))
        x = np.linspace(log_mass[0], log_mass[-1], layers)
        for field in ("temp", "xne", "xna", "rho"):
            model[field] = np.exp(PchipInterpolator(log_mass, np.log(base[field]))(x))
        model.rhox = np.exp(x)
        model.geom = "SPH"
        model.radius = 1e11
        model.height = 2e9 * (1 - (x-log_mass[0])/(log_mass[-1]-log_mass[0])) + 1e7
        model.wlstd = 5000.0
        model.depth = "RHOX"
        return model

    def prepare(model, scattering):
        dll.InputLineList(lines)
        dll.InputModel(5750.0, 2.5, 2.0, model)
        dll.InputAbund(Abund(0.0, "asplund2009"))
        # Preserve the supplied densities while calculating EOS populations.
        dll.Ionization(0)
        dll.SetContinuumScatteringSourceMode(scattering)
        dll.InputWaveRange(wave[0]-0.01, wave[-1]+0.01)
        dll.Opacity()

    def spectrum(model, scattering):
        prepare(model, scattering)
        _, _, intensity, continuum = dll.Transf(mu, wave=wave)
        if not np.all(np.isfinite(intensity)) or not np.all(continuum > 0):
            raise AssertionError("Non-finite intensity or non-positive continuum")
        return intensity, continuum

    try:
        # Compute the transformation on an independent dense grid, then use
        # the same mapping for every transfer resolution. Fine integration
        # in log(m) reduces quadrature error relative to coarse-grid DELO.
        dense = atmosphere(288)
        prepare(dense, False)
        _, _, chi = dll.GetContinuumOpacityComponents(5000.0)
        x = np.log(dense.rhox)
        fine_x = np.linspace(x[0], x[-1], 30001)
        fine_chi = np.exp(PchipInterpolator(x, np.log(chi))(fine_x))
        tau0 = chi[0] * dense.rhox[0]
        fine_tau = tau0 + cumulative_trapezoid(
            fine_chi * np.exp(fine_x), x=fine_x, initial=0,
        )
        tau_from_log_mass = PchipInterpolator(fine_x, fine_tau)
        rows = []
        for scattering in (False, True):
            for layers in (24, 48, 96, 192, 288):
                model = atmosphere(layers)
                model.tau = tau_from_log_mass(np.log(model.rhox))
                ir, cr = spectrum(model, scattering)
                model.depth = "TAU"
                it, ct = spectrum(model, scattering)
                # Near-zero continuum on a grazing ray is a poor denominator
                # for line differences. Scale by that ray's peak intensity.
                errors = np.max(np.abs(ir-it), axis=1) / np.max(np.abs(ir), axis=1)
                continuum_error = np.max(np.abs(cr-ct)/cr, axis=1)
                row = {
                    "scattering_source": int(scattering), "layers": layers,
                    "peak_scaled_intensity_error_by_mu": errors.tolist(),
                    "max_peak_scaled_intensity_error": float(np.max(errors)),
                    "relative_continuum_error_by_mu": continuum_error.tolist(),
                    "max_relative_continuum_error": float(np.max(continuum_error)),
                }
                rows.append(row)
                np.savez_compressed(
                    args.out / f"scattering-{int(scattering)}-layers-{layers}.npz",
                    wave=wave, mu=mu, intensity_rhox=ir, intensity_tau=it,
                    continuum_rhox=cr, continuum_tau=ct,
                    rhox=model.rhox, tau=model.tau, height=model.height,
                )
                print(json.dumps(row), flush=True)
        for scattering in (0, 1):
            cases = [row for row in rows if row["scattering_source"] == scattering]
            assert cases[-1]["max_peak_scaled_intensity_error"] < 1e-3
            assert cases[-1]["max_relative_continuum_error"] < 2e-3
            assert cases[-1]["max_peak_scaled_intensity_error"] < cases[0]["max_peak_scaled_intensity_error"] / 20
        report = {
            "native_version": dll.SMELibraryVersion(),
            "native_library": str(Path(get_full_libfile()).resolve()),
            "native_sha256": hashlib.sha256(Path(get_full_libfile()).read_bytes()).hexdigest(),
            "atmosphere_sha256": hashlib.sha256(args.atmosphere.read_bytes()).hexdigest(),
            "reference_wavelength_A": 5000.0, "reference_layers": 288,
            "quadrature_points": 30001, "ionization_flags": 0,
            "teff_K": 5750.0, "logg": 2.5, "monh": 0.0,
            "radius_cm": 1e11, "mu": mu.tolist(),
            "wave_A": [float(wave[0]), float(wave[-1]), len(wave)],
            "rows": rows,
        }
        (args.out / "report.json").write_text(json.dumps(report, indent=2) + "\n")
    finally:
        dll.SetContinuumScatteringSourceMode(False)
        dll.ResetDepartureCoefficients()


if __name__ == "__main__":
    main()
