# -*- coding: utf-8 -*-
"""Validate spherical interpolation against an off-grid MARCS reference model."""
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from pysme.atmosphere.interpolation import AtmosphereInterpolator

from .test_largefilestorage import lfs_atmo, skipif_lfs

SPHERICAL_GRID_NAME = "marcs2014.sav"
# Exact off-grid MARCS model computed by Nils Ryde using the MARCS code.
REFERENCE_LIS_PATH = Path(__file__).parent / "5035g2.96z-0.74m1.0t1_new.lis"
TARGET_TEFF = 5035.0
TARGET_LOGG = 2.96
TARGET_MONH = -0.74


def _find_marker_lines(path):
    """Locate the 0-indexed line numbers MARCS' fixed-format .lis tables start at."""
    lineone = linetwo = lineradius = None
    with open(path) as f:
        for i, line in enumerate(f):
            if "Spherical model with  Radius" in line:
                lineradius = i
            elif "M o d e l   A t m o s p h e r e     (cgs units)" in line:
                lineone = i + 4
            elif (
                "T h e r m o d y n a m i c a l   q u a n t i t i e s  and   "
                "C o n v e c t i o n   (cgs units)" in line
            ):
                linetwo = i + 4
    return lineone, linetwo, lineradius


def _read_marcs_lis_reference(path, ndep=56):
    """Parse a MARCS ``.lis`` model file into the quantities used by PySME."""
    lineone, linetwo, lineradius = _find_marker_lines(path)

    df1 = pd.read_csv(path, skiprows=lineone, header=None, nrows=ndep, sep=r"\s+")
    df1.columns = [
        "K", "TauRoss", "Tau5000", "GeomDepth", "T", "Pe", "Pg", "Prad",
        "Pturb", "KappaRoss", "K2",
    ]

    df2 = pd.read_csv(path, skiprows=linetwo, header=None, nrows=ndep, sep=r"\s+")
    df2.columns = [
        "K", "TauRoss_2", "Density", "Mu", "Cp", "Cv", "AdGrad", "Q",
        "SoundVel", "ConvVel", "FconvF", "K2_2",
    ]
    merged = df1.merge(df2, on="K")

    radiustemp = pd.read_csv(path, skiprows=lineradius, header=None, nrows=1, sep=r"or")
    radiustemp.columns = ["temp1", "radiusactual"]
    radius = float(radiustemp["radiusactual"].values[0].replace("cm", ""))

    kb = 1.380649e-16
    temp = merged["T"].to_numpy(dtype=float)
    tau = merged["Tau5000"].to_numpy(dtype=float)
    tauross = merged["TauRoss"].to_numpy(dtype=float)
    rho = merged["Density"].to_numpy(dtype=float)
    xne = merged["Pe"].to_numpy(dtype=float) / (kb * temp)
    xna = merged["Pg"].to_numpy(dtype=float) / (kb * temp) - xne
    height = -merged["GeomDepth"].to_numpy(dtype=float)
    kappaross = merged["KappaRoss"].to_numpy(dtype=float)

    # Mass column density, by the same trapezoidal d(tau)/kappa integration
    # used to derive it for this reference model outside of PySME.
    rhox = np.zeros_like(tau)
    for i in range(1, len(tauross)):
        delta_tau = tauross[i] - tauross[i - 1]
        avg_inv_kappa = 0.5 * (1.0 / kappaross[i] + 1.0 / kappaross[i - 1])
        rhox[i] = rhox[i - 1] + delta_tau * avg_inv_kappa

    return dict(
        temp=temp, tau=tau, rho=rho, xne=xne, xna=xna, height=height,
        radius=radius, rhox=rhox,
    )


@skipif_lfs
def test_spherical_interp_matches_offgrid_marcs_reference(lfs_atmo):
    interpolator = AtmosphereInterpolator(
        depth="RHOX", interp="RHOX", geom="SPH", lfs_atmo=lfs_atmo
    )
    atmo = interpolator.interp_atmo_grid(
        SPHERICAL_GRID_NAME, TARGET_TEFF, TARGET_LOGG, TARGET_MONH
    )
    reference = _read_marcs_lis_reference(REFERENCE_LIS_PATH)

    # Compare on the reference's rhox range only; interp_atmo_grid's own
    # depth scale need not span exactly the same rhox extent.
    in_range = (atmo.rhox >= reference["rhox"].min()) & (
        atmo.rhox <= reference["rhox"].max()
    )
    assert np.count_nonzero(in_range) > 0

    # An off-grid MARCS model isn't expected to match the interpolated grid
    # point-for-point (this reference was computed with slightly different
    # physics, e.g. CN-cycling, than marcs2014.sav), so this checks overall
    # closeness (RMS error) rather than a tight per-point tolerance. The
    # thresholds are calibrated with margin above what combined height+radius
    # handling actually measures against this reference (RMS ~4.3e7 cm,
    # max ~1.8e8 cm) - both well under half of what separate height/radius
    # handling measures against the same reference (RMS ~2.2e8 cm,
    # max ~3.3e8 cm), confirming the combined approach is a real improvement.
    expected_height = np.interp(atmo.rhox[in_range], reference["rhox"], reference["height"])
    height_diff = atmo.height[in_range] - expected_height
    assert np.sqrt(np.mean(height_diff ** 2)) < 1e8
    assert np.max(np.abs(height_diff)) < 3e8

    expected_temp = np.interp(atmo.rhox[in_range], reference["rhox"], reference["temp"])
    assert np.allclose(atmo.temp[in_range], expected_temp, rtol=0.05)
