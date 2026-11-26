# -*- coding: utf-8 -*-
"""Regression test for the Stark profiles of Brackett lines (lines without Stehle
tables, computed with HLINOP).

Since SMElib 6.13.22 the broadening mechanisms are convolved for all such lines
by default (h_line_profile="convolved"); h_line_profile="legacy" restores the
earlier additive profiles, whose Stark areas were 1.7-2.1.  The earlier opt-in
switch h_stark_convolution (Brackett m >= 10 only) no longer has an effect.

Each case runs in an isolated subprocess so native-library state cannot leak
between cases.  The thresholds guard against silent regressions:
  - legacy reproduces the earlier fixed reference EW
  - the default (convolved) reproduces its reference EW and differs from legacy
  - h_stark_convolution does not change the result
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
    # name: (center_air, half_window, reference_legacy_EW, reference_convolved_EW)
    # legacy: the earlier fixed reference (unchanged); convolved: SMElib 6.13.22
    "br11": (16806.538, 5.0, 0.92623631, 0.44376523),
    "br10": (17362.112, 5.0, 1.21278839, 0.65029839),
}

# The convolved profiles must be clearly different from legacy.
MIN_CONV_VS_LEGACY_DIFF = 0.05  # EW difference in Angstrom
CONV_EW_REL_TOL = 0.005  # 0.5%
LEGACY_EW_REL_TOL = 0.005  # 0.5%

WORKER = r'''
import os, sys
import numpy as np
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
profile, stark = {profile!r}, {stark!r}
if profile is not None:
    sme.h_line_profile = profile
if stark is not None:
    sme.h_stark_convolution = stark
out = synthesize_spectrum(sme, linelist_mode="all")
flux = np.asarray(out.synth[0], dtype=float)
ew = float(np.trapezoid(1.0 - flux, wave))
print("{{:.8f}}".format(ew))
'''


def _synthesize_ew(profile, center: float, half: float, stark=None) -> float:
    # None leaves the SME_Structure defaults untouched.
    script = WORKER.format(
        profile=profile,
        stark=stark,
        sun=SUN,
        linelist=str(LINELIST),
        center=center,
        half=half,
    )
    env = os.environ.copy()
    env.pop("PYSME_H_STARK_CONVOLUTION", None)  # no longer read by SMElib
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
def test_brackett_profiles(line):
    center, half, ref_legacy_ew, ref_conv_ew = LINES[line]

    ew_legacy = _synthesize_ew("legacy", center, half)
    ew_conv = _synthesize_ew(None, center, half)   # default: convolved

    legacy_rel = abs(ew_legacy - ref_legacy_ew) / ref_legacy_ew
    assert legacy_rel < LEGACY_EW_REL_TOL, (
        f"{line}: legacy EW {ew_legacy:.8f} deviates {legacy_rel*100:.2f}% "
        f"from fixed reference {ref_legacy_ew:.8f}"
    )
    assert abs(ew_conv - ew_legacy) > MIN_CONV_VS_LEGACY_DIFF, (
        f"{line}: convolved ({ew_conv:.4f}) should differ from legacy "
        f"({ew_legacy:.4f})"
    )
    rel = abs(ew_conv - ref_conv_ew) / ref_conv_ew
    assert rel < CONV_EW_REL_TOL, (
        f"{line}: convolved EW {ew_conv:.4f} deviates {rel*100:.1f}% "
        f"from reference {ref_conv_ew:.4f}"
    )


def test_h_stark_convolution_has_no_effect():
    """The deprecated switch h_stark_convolution no longer changes Br11."""
    center, half, _, _ = LINES["br11"]
    ew_default = _synthesize_ew(None, center, half)
    ew_switch = _synthesize_ew(None, center, half, stark="convolution")
    assert ew_switch == ew_default
