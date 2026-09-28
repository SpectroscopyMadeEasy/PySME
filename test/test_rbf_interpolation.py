# -*- coding: utf-8 -*-
"""Tests for the RBF-based atmosphere-grid interpolation path (RbfGrid/interpolate_RBF)."""
import numpy as np
import pytest

from pysme.abund import Abund, elements_dict
from pysme.atmosphere.atmosphere import AtmosphereError, AtmosphereGrid
from pysme.atmosphere.interpolation import AtmosphereInterpolator
from pysme.atmosphere.savfile import SavFile

from .test_largefilestorage import lfs_atmo, skipif_lfs

# Approximate solar photospheric abundances (H=12 scale), in the standard
# element order (H, He, Li, Be, B, C, N, O, ...). Illustrative values in the
# same ballpark as Grevesse & Sauval (1998)/Asplund et al. (2009) - not a
# citation-grade table, just realistic enough that test fixtures don't look
# like an arbitrary constant. NaN marks elements with no well-established
# solar value (unstable isotopes Tc/Pm, or heavy elements beyond Bi besides
# Th/U), matching how a real solar pattern leaves them untracked.
_SOLAR_H12 = np.array([
    12.00, 10.93,                                                  # H, He
    1.05, 1.38, 2.70, 8.43, 7.83, 8.69, 4.56, 7.93,                 # Li-Ne
    6.24, 7.60, 6.45, 7.51, 5.41, 7.12, 5.50, 6.40,                 # Na-Ar
    5.03, 6.34, 3.15, 4.95, 3.93, 5.64, 5.43, 7.50,                 # K-Fe
    4.99, 6.22, 4.19, 4.56, 3.04, 3.65, 2.30, 3.34,                 # Co-Se
    2.54, 3.25, 2.52, 2.87, 2.21, 2.58, 1.46, 1.88,                 # Br-Mo
    np.nan, 1.75, 0.91, 1.57, 0.94, 1.71, 0.80, 2.04,               # Tc-Sn
    1.01, 2.18, 1.55, 2.24, 1.08, 2.18, 1.10, 1.58,                 # Sb-Ce
    0.72, 1.42, np.nan, 0.96, 0.52, 1.07, 0.30, 1.10,               # Pr-Dy
    0.48, 0.92, 0.10, 0.84, 0.10, 0.85, -0.12, 0.85,                # Ho-W
    0.26, 1.40, 1.38, 1.62, 0.92, 1.17, 0.90, 1.75,                 # Re-Pb
    0.65, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, 0.02,     # Bi-Th
    np.nan, -0.54, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan,  # Pa-Cf
    np.nan,                                                         # Es
])
assert len(_SOLAR_H12) == 99


def _make_synthetic_grid(
    source, teffs, loggs, monhs, temp_value, ndep=3, geom="PP",
    opflag_value=1, wlstd_value=5000.0, abund_shift=0.0,
):
    """
    Build a tiny in-memory AtmosphereGrid spanning the Cartesian product of
    teffs x loggs x monhs, with every depth point of ``temp`` set to
    ``temp_value`` and every other field set to a fixed, physically-harmless
    placeholder - enough for RbfGrid to build an interpolator, using the same
    field-assignment convention SavFile uses when loading a real grid.

    The abundance pattern is realistic solar (``_SOLAR_H12``) with
    ``abund_shift`` added to the metals (index 2+) - a constant, or, for
    tests that need abundance to vary with metallicity, a callable of monh
    (mimicking a realized, metallicity-baked-in pattern: solar + monh).
    Stored in "sme" format (via ``Abund.totype``) like a real grid, so
    reading it back through ``_realized_abund_h12`` round-trips correctly.
    """
    combos = [(t, g, m) for t in teffs for g in loggs for m in monhs]
    natmo = len(combos)
    grid = AtmosphereGrid(natmo, ndep)
    grid.source = source
    grid.geom = geom
    grid["teff"] = [c[0] for c in combos]
    grid["logg"] = [c[1] for c in combos]
    grid["monh"] = [c[2] for c in combos]
    grid["vturb"] = 0.0
    grid["lonh"] = 0.0
    grid["wlstd"] = wlstd_value
    grid["radius"][:] = 1 if geom == "PP" else 1e11
    grid["temp"] = temp_value
    grid["tau"] = np.tile(np.linspace(1e-4, 1e-1, ndep), (natmo, 1))
    grid["rhox"] = np.tile(np.linspace(1e-3, 1.0, ndep), (natmo, 1))
    grid["rho"] = 1.0
    grid["xna"] = 1.0
    grid["xne"] = 1.0
    grid["opflag"] = opflag_value
    for i, (_, _, m) in enumerate(combos):
        shift = abund_shift(m) if callable(abund_shift) else abund_shift
        intended_h12 = _SOLAR_H12.copy()
        intended_h12[2:] += shift
        grid["abund"][i] = Abund.totype(intended_h12, "sme", raw=True)
    if geom == "SPH":
        grid["height"] = np.tile(np.linspace(0.0, 1.0, ndep), (natmo, 1))
    return grid


def test_rbf_cache_invalidated_on_source_change():
    teffs, loggs, monhs = [4900.0, 5100.0], [3.8, 4.2], [-0.2, 0.2]
    grid_a = _make_synthetic_grid("grid_a", teffs, loggs, monhs, temp_value=5000.0)
    grid_b = _make_synthetic_grid("grid_b", teffs, loggs, monhs, temp_value=6000.0)

    interpolator = AtmosphereInterpolator(interp="RBF")
    atmo_a = interpolator.interp_atmo_grid(grid_a, 5000.0, 4.0, 0.0)
    # Same interpolator instance, different grid: must not reuse grid_a's cached RbfGrid.
    atmo_b = interpolator.interp_atmo_grid(grid_b, 5000.0, 4.0, 0.0)

    fresh_interpolator = AtmosphereInterpolator(interp="RBF")
    atmo_b_fresh = fresh_interpolator.interp_atmo_grid(grid_b, 5000.0, 4.0, 0.0)

    assert np.allclose(atmo_a.temp, 5000.0)
    assert np.allclose(atmo_b.temp, atmo_b_fresh.temp)
    assert not np.allclose(atmo_b.temp, atmo_a.temp)


def test_rbf_handles_single_valued_metallicity_axis():
    # A single metallicity value (like the real spherical MARCS grid) must
    # not crash the per-axis step-size computation in initialize_gridpoints.
    grid = _make_synthetic_grid(
        "single_monh_grid", teffs=[4900.0, 5100.0], loggs=[3.8, 4.2], monhs=[0.0],
        temp_value=5000.0,
    )
    interpolator = AtmosphereInterpolator(interp="RBF")
    atmo = interpolator.interp_atmo_grid(grid, 5000.0, 4.0, 0.0)
    assert np.all(np.isfinite(atmo.temp))


def test_rbf_preserves_opflag_and_wlstd():
    # opflag and wlstd deliberately differ here from Atmo()'s hardcoded
    # defaults ([1]*20, 5000.0).
    teffs, loggs, monhs = [4900.0, 5100.0], [3.8, 4.2], [-0.2, 0.2]
    grid = _make_synthetic_grid(
        "metadata_grid", teffs, loggs, monhs, temp_value=5000.0,
        opflag_value=0, wlstd_value=6000.0,
    )
    interpolator = AtmosphereInterpolator(interp="RBF")
    atmo = interpolator.interp_atmo_grid(grid, 5000.0, 4.0, 0.0)

    assert np.array_equal(atmo.opflag, np.zeros(20, dtype=atmo.opflag.dtype))
    assert atmo.wlstd == 6000.0


def test_rbf_abundance_stays_consistent_with_interpolated_metallicity():
    # Every grid point's abundance is solar shifted by exactly its own monh
    # (a stand-in for a realized, metallicity-baked-in pattern, as confirmed
    # against the real marcs2014.sav grid). teff/logg are irrelevant to the
    # abundance here, only monh matters.
    teffs, loggs, monhs = [4900.0, 5100.0], [3.8, 4.2], [-2.0, -1.0, 0.0]
    grid = _make_synthetic_grid(
        "abundance_grid", teffs, loggs, monhs, temp_value=5000.0,
        abund_shift=lambda monh: monh,
    )
    interpolator = AtmosphereInterpolator(interp="RBF")

    # Off-grid metallicity, roughly midway between two grid values.
    query_monh = -1.4
    atmo = interpolator.interp_atmo_grid(grid, 5000.0, 4.0, query_monh)

    # atmo.monh itself is only approximately the query (RBF is a fit, not
    # exact interpolation at off-grid points), so compare the abundance
    # against a pattern built the same way at atmo's *own* interpolated
    # monh - this is the actual bug fix under test: the returned abundance
    # must reflect the atmosphere's own reported metallicity, not whichever
    # grid point happened to be nearest.
    expected_fe = _SOLAR_H12[elements_dict["Fe"]] + atmo.monh

    assert np.isclose(atmo.abund.get_pattern_abundance("Fe"), expected_fe, atol=1e-6)


def test_rbf_interpolates_spherical_geometry():
    # Every other RBF test uses the default geom="PP", so none of them
    # exercise the SPH-only vtags (height combined with radius, per
    # initialize_tags/to_interp_space_vector). Use a non-degenerate 2x2x2
    # grid so this isolates the SPH-specific readout path from the
    # single-value-axis handling covered separately above.
    teffs, loggs, monhs = [4500.0, 4700.0], [2.0, 2.4], [-0.5, 0.0]
    ndep = 4
    grid = _make_synthetic_grid(
        "spherical_grid", teffs, loggs, monhs, temp_value=4600.0,
        ndep=ndep, geom="SPH",
    )
    interpolator = AtmosphereInterpolator(interp="RBF", geom="SPH")
    atmo = interpolator.interp_atmo_grid(grid, 4600.0, 2.2, -0.25)

    assert np.all(np.isfinite(atmo.height))
    assert len(atmo.height) == ndep
    assert atmo.radius > 0
    # The helper uses radius=1e11 for SPH grids (vs. radius=1 for PP) - a
    # sanity check that the interpolated radius lands in that same regime,
    # rather than e.g. silently falling back to a PP-shaped result.
    assert np.isclose(atmo.radius, 1e11, rtol=0.5)


def test_rbf_out_of_domain_raises_with_error_policy():
    teffs, loggs, monhs = [4900.0, 5100.0], [3.8, 4.2], [-0.2, 0.2]
    grid = _make_synthetic_grid("bounds_grid", teffs, loggs, monhs, temp_value=5000.0)
    interpolator = AtmosphereInterpolator(interp="RBF")

    with pytest.raises(AtmosphereError):
        interpolator.interp_atmo_grid(
            grid, 50000.0, 4.0, 0.0, interpolation_policy="error"
        )


def test_rbf_out_of_domain_extrapolates_with_allow_policy():
    teffs, loggs, monhs = [4900.0, 5100.0], [3.8, 4.2], [-0.2, 0.2]
    grid = _make_synthetic_grid("bounds_grid", teffs, loggs, monhs, temp_value=5000.0)
    interpolator = AtmosphereInterpolator(interp="RBF")

    # Default policy ("allow"): still returns a (extrapolated) result rather than raising.
    atmo = interpolator.interp_atmo_grid(grid, 50000.0, 4.0, 0.0)
    assert np.all(np.isfinite(atmo.temp))


def _real_depth_count(atmo):
    """
    marcs2014.sav pads every atmosphere's per-depth arrays to a fixed max
    length (72) with trailing zeros/NaNs beyond that atmosphere's real depth
    count (56, for every point checked so far). RbfGrid.initialize_depth
    already strips this globally (from the first grid atmosphere) before
    fitting; grid-fetched comparison values need the same stripping applied
    per-atmosphere before comparing against RBF's (already-truncated) output.
    """
    tau = np.asarray(atmo.tau)
    if np.isnan(tau[-1]):
        return int(np.sum(~np.isnan(tau)))
    if tau[-1] == 0:
        return int(np.sum(tau > 0.0))
    return len(tau)


@skipif_lfs
def test_rbf_interpolation_matches_real_grid(lfs_atmo):
    """
    Tie-together sanity check against a real MARCS grid (marcs2014.sav),
    covering both the atmosphere-structure interpolation and the abundance
    interpolation together. Chosen over the smaller, mono-metallic
    marcs2012s_t1.0.sav for this test because it's already used elsewhere in
    this project (no extra download) and lets abundance-interpolation
    smoothness be checked at the same time.

    Note: the *absolute* abundance values from this grid are known to be
    wrong (a separate, pre-existing Abund/Atmosphere double-counting bug,
    unrelated to RBF/PR #22 - see the write-up sent to the maintainer). This
    test only checks that RBF reproduces the (still-realized) pattern
    exactly at a grid point and interpolates it smoothly off-grid, not that
    the values are physically correct.
    """
    grid = SavFile(lfs_atmo.get("marcs2014.sav"), source="marcs2014.sav", lfs=lfs_atmo)
    interpolator = AtmosphereInterpolator(interp="RBF", geom="SPH", lfs_atmo=lfs_atmo)
    logg, monh = 2.5, -1.0

    # Exact grid point: RBFInterpolator (smoothing=0) is an exact interpolant
    # at its own input points (the query is its own nearest neighbor), so
    # this should reproduce the stored model almost exactly - the real-data
    # analog of test_atmospheres.py::test_grid_point, without that test's
    # depth-index shift (RBF doesn't clip a top point the way the TAU/
    # pairwise path does).
    teff = 5000.0
    atmo_interp = interpolator.interp_atmo_grid(grid, teff, logg, monh)
    atmo_grid = grid.get(teff, logg, monh)
    n = _real_depth_count(atmo_grid)

    assert np.allclose(atmo_interp.temp, np.asarray(atmo_grid.temp)[:n], rtol=1e-6)
    assert np.allclose(
        atmo_interp.height, np.asarray(atmo_grid.height)[:n], rtol=1e-6, atol=1e-3
    )
    assert np.isclose(atmo_interp.radius, atmo_grid.radius, rtol=1e-6)
    assert np.isclose(
        atmo_interp.abund.get_pattern_abundance("Si"),
        atmo_grid.abund.get_pattern_abundance("Si"),
        rtol=1e-6,
    )

    # Off-grid teff, bracketed by two real grid points at the same logg/monh.
    # No independent reference exists for an arbitrary off-grid point, so
    # this stays a plausibility check for both the structure and the
    # abundance: finite, still hotter with depth, and landing close to (not
    # necessarily exactly inside) the bracketing points' own values.
    lower = grid.get(4750.0, logg, monh)
    upper = grid.get(5000.0, logg, monh)
    atmo_between = interpolator.interp_atmo_grid(grid, 4875.0, logg, monh)

    n_lo, n_up = _real_depth_count(lower), _real_depth_count(upper)
    assert n_lo == n_up == len(atmo_between.temp)
    lo_temp = np.asarray(lower.temp)[:n_lo]
    up_temp = np.asarray(upper.temp)[:n_up]

    assert np.all(np.isfinite(atmo_between.temp))
    assert np.all(np.isfinite(atmo_between.height))
    assert np.all(np.diff(atmo_between.temp) >= -1e-6 * np.max(atmo_between.temp))

    lo_bound = np.minimum(lo_temp, up_temp)
    hi_bound = np.maximum(lo_temp, up_temp)
    margin = 0.1 * (hi_bound - lo_bound)
    assert np.all(atmo_between.temp >= lo_bound - margin)
    assert np.all(atmo_between.temp <= hi_bound + margin)

    si_between = atmo_between.abund.get_pattern_abundance("Si")
    si_lo = lower.abund.get_pattern_abundance("Si")
    si_up = upper.abund.get_pattern_abundance("Si")
    assert np.isfinite(si_between)
    lo_si, hi_si = sorted((si_lo, si_up))
    si_margin = max(0.1 * (hi_si - lo_si), 0.05)
    assert lo_si - si_margin <= si_between <= hi_si + si_margin


@skipif_lfs
def test_rbf_interpolation_matches_real_grid_pp(lfs_atmo):
    """
    Real-grid tie-together test for the plane-parallel (PP) counterpart of
    test_rbf_interpolation_matches_real_grid. Uses marcs2012p_t1.0.sav
    (already used by test_atmospheres.py::test_grid_point, so no new
    download) rather than reusing marcs2014.sav's approach: it has a
    uniform 56-point depth array with zero padding (confirmed via a
    full-grid scan), so no truncation helper is needed here, and it still
    has real Si variation across its 15 metallicities for a meaningful
    abundance-smoothness check. PP has no height/radius vtag, so those
    checks from the SPH test don't apply here.
    """
    grid = SavFile(
        lfs_atmo.get("marcs2012p_t1.0.sav"), source="marcs2012p_t1.0.sav", lfs=lfs_atmo
    )
    interpolator = AtmosphereInterpolator(interp="RBF", geom="PP", lfs_atmo=lfs_atmo)
    logg, monh = 4.0, -1.0

    # Exact grid point: RBFInterpolator (smoothing=0) is an exact interpolant
    # at its own input points, so this should reproduce the stored model
    # (structure and abundance) almost exactly.
    teff = 5000.0
    atmo_interp = interpolator.interp_atmo_grid(grid, teff, logg, monh)
    atmo_grid = grid.get(teff, logg, monh)

    assert np.allclose(atmo_interp.temp, atmo_grid.temp, rtol=1e-6)
    assert np.isclose(
        atmo_interp.abund.get_pattern_abundance("Si"),
        atmo_grid.abund.get_pattern_abundance("Si"),
        rtol=1e-6,
    )

    # Off-grid teff, bracketed by two real grid points at the same logg/monh.
    # No independent reference exists for an arbitrary off-grid point, so
    # this stays a plausibility check: finite, still hotter with depth, and
    # landing close to (not necessarily exactly inside) the bracketing
    # points' own values.
    lower = grid.get(4750.0, logg, monh)
    upper = grid.get(5000.0, logg, monh)
    atmo_between = interpolator.interp_atmo_grid(grid, 4875.0, logg, monh)

    assert np.all(np.isfinite(atmo_between.temp))
    assert np.all(np.diff(atmo_between.temp) >= -1e-6 * np.max(atmo_between.temp))

    lo_bound = np.minimum(lower.temp, upper.temp)
    hi_bound = np.maximum(lower.temp, upper.temp)
    margin = 0.1 * (hi_bound - lo_bound)
    assert np.all(atmo_between.temp >= lo_bound - margin)
    assert np.all(atmo_between.temp <= hi_bound + margin)

    si_between = atmo_between.abund.get_pattern_abundance("Si")
    si_lo = lower.abund.get_pattern_abundance("Si")
    si_up = upper.abund.get_pattern_abundance("Si")
    assert np.isfinite(si_between)
    lo_si, hi_si = sorted((si_lo, si_up))
    si_margin = max(0.1 * (hi_si - lo_si), 0.05)
    assert lo_si - si_margin <= si_between <= hi_si + si_margin
