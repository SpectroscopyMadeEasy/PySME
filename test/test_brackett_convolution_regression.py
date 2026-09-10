# -*- coding: utf-8 -*-
"""Regression test for the opt-in Brackett shared-support Stark convolution.

The convolution is OFF by default (legacy additive construction).  Setting
PYSME_H_STARK_CONVOLUTION=convolution enables the shared-support convolution
profile for Brackett lines with m>=10 (Br10 and above).

Each mode runs in an isolated subprocess so process-global environment and
native-library state cannot leak between the legacy and convolution cases.

Reference EW values were validated against the analysis prototype and Korg
(see analysis/h_occ_validation/).  The thresholds here only guard against
silent regressions:
  - convolution must differ from legacy (the whole point of the switch)
  - convolution must reproduce the accepted reference EW within tolerance
"""
from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
LINELIST = Path(__file__).with_name("brackett_br10_br11.lin")

SUN = {"teff": 5772.0, "logg": 4.44, "monh": 0.0, "vmic": 1.0, "vmac": 2.5, "vsini": 1.8}

LINES = {
    # name: (center_air, half_window, reference_legacy_EW, reference_conv_EW)
    "br11": (16806.538, 5.0, 0.92623631, 0.35522724),
    "br10": (17362.112, 5.0, 1.21278839, 0.56173425),
}

# Convolution must be clearly different from legacy (well above numerical noise).
MIN_CONV_VS_LEGACY_DIFF = 0.05  # EW difference in Angstrom
# Convolution EW must reproduce the accepted reference within this tolerance.
CONV_EW_REL_TOL = 0.005  # 0.5%
LEGACY_EW_REL_TOL = 0.005  # 0.5%

WORKER = r'''
import os, sys
from pathlib import Path
import numpy as np
sys.meta_path[:] = [
    finder for finder in sys.meta_path
    if finder.__class__.__module__ != "_editable_skbc_pysme_astro"
]
sys.path.insert(0, str(Path({root!r}) / "src"))
from pysme.abund import Abund
from pysme.linelist.vald import ValdFile
from pysme.sme import SME_Structure
from pysme.synthesize import synthesize_spectrum
sme = SME_Structure()
for k, v in {sun!r}.items():
    setattr(sme, k, v)
sme.iptype = "gauss"
sme.ipres = 300000.0
sme.atmo.method = "grid"
sme.atmo.source = "marcs2012.sav"
sme.abund = Abund(pattern="asplund2009", monh=0.0)
sme.linelist = ValdFile({linelist!r})
wave = np.arange({center} - {half}, {center} + {half} + 0.01, 0.02)
sme.wave = [wave]
sme.normalize_by_continuum = True
sme.vrad_flag = "none"
sme.cscale_flag = "none"
sme.normalize_resample_mode = "ratio"
mode = {mode!r}
if mode is not None:
    sme.h_stark_convolution = mode
out = synthesize_spectrum(sme, linelist_mode="all")
flux = np.asarray(out.synth[0], dtype=float)
ew = float(np.trapezoid(1.0 - flux, wave))
print("{{:.8f}}".format(ew))
'''


def _synthesize_ew(mode: str, center: float, half: float) -> float:
    # Explicit modes use the public API field.  None leaves the SME_Structure
    # default untouched for the default-path test.
    script = WORKER.format(
        root=str(ROOT),
        mode=mode,
        sun=SUN,
        linelist=str(LINELIST),
        center=center,
        half=half,
    )
    env = os.environ.copy()
    env.pop("PYSME_H_STARK_CONVOLUTION", None)  # rely on the API field, not env
    env["PYSME_RESAMPLE_NORM_MODE"] = "ratio"
    proc = subprocess.run(
        [sys.executable, "-c", script],
        cwd=str(ROOT),
        env=env,
        capture_output=True,
        text=True,
    )
    if proc.returncode != 0:
        pytest.fail(
            f"Brackett convolution synthesis failed: {proc.stderr[-1000:]}",
            pytrace=False,
        )
    return float(proc.stdout.strip().splitlines()[-1])


@pytest.mark.parametrize("line", ["br11", "br10"])
def test_brackett_convolution_switch(line):
    center, half, ref_legacy_ew, ref_conv_ew = LINES[line]

    ew_legacy = _synthesize_ew("legacy", center, half)
    ew_conv = _synthesize_ew("convolution", center, half)

    legacy_rel = abs(ew_legacy - ref_legacy_ew) / ref_legacy_ew
    assert legacy_rel < LEGACY_EW_REL_TOL, (
        f"{line}: legacy EW {ew_legacy:.8f} deviates {legacy_rel*100:.2f}% "
        f"from fixed reference {ref_legacy_ew:.8f}"
    )

    # The convolution must actually change the profile vs legacy.
    assert abs(ew_conv - ew_legacy) > MIN_CONV_VS_LEGACY_DIFF, (
        f"{line}: convolution ({ew_conv:.4f}) should differ from legacy "
        f"({ew_legacy:.4f})"
    )

    # Convolution EW must match the accepted reference.
    rel = abs(ew_conv - ref_conv_ew) / ref_conv_ew
    assert rel < CONV_EW_REL_TOL, (
        f"{line}: convolution EW {ew_conv:.4f} deviates {rel*100:.1f}% "
        f"from reference {ref_conv_ew:.4f}"
    )


def test_brackett_convolution_default_is_legacy():
    """Without the API field set (None), Br11 must match the legacy profile."""
    center, half, _, _ = LINES["br11"]

    # mode=None leaves SME_Structure at its explicit PySME default (legacy).
    ew_default = _synthesize_ew(None, center, half)

    ew_legacy = _synthesize_ew("legacy", center, half)
    assert abs(ew_default - ew_legacy) < 1e-6, (
        f"default EW {ew_default:.6f} must equal legacy EW {ew_legacy:.6f}"
    )
