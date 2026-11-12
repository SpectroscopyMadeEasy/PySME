# -*- coding: utf-8 -*-
"""Tests for RBF atmosphere-grid interpolation (RbfGrid, interpolate_RBF)."""
import warnings

import numpy as np
import pytest

from pysme.abund import Abund, elements_dict
from pysme.atmosphere.atmosphere import AtmosphereError, AtmosphereGrid
from pysme.atmosphere.interpolation import AtmosphereInterpolator, RbfGrid
from pysme.atmosphere.savfile import SavFile

from .test_largefilestorage import lfs_atmo, skipif_lfs

# Approximate solar abundances (H=12), roughly Grevesse & Sauval (1998) /
# Asplund et al. (2009); realistic test values, not a reference table. NaN marks
# elements without a solar value, as in a real solar pattern.
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
    Small in-memory AtmosphereGrid over teffs x loggs x monhs, with temp set to
    temp_value and placeholder values for the other fields.

    The abundances are _SOLAR_H12 with abund_shift added to the metals:
    a constant, or a function of monh (e.g. solar + monh, like a realized
    pattern). They are stored in "sme" format, like a real grid.
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
    # Same interpolator, different grid: must not reuse grid_a's RbfGrid
    atmo_b = interpolator.interp_atmo_grid(grid_b, 5000.0, 4.0, 0.0)

    fresh_interpolator = AtmosphereInterpolator(interp="RBF")
    atmo_b_fresh = fresh_interpolator.interp_atmo_grid(grid_b, 5000.0, 4.0, 0.0)

    assert np.allclose(atmo_a.temp, 5000.0)
    assert np.allclose(atmo_b.temp, atmo_b_fresh.temp)
    assert not np.allclose(atmo_b.temp, atmo_a.temp)


def _single_monh_grid():
    # One metallicity, like the spherical MARCS grid, with a varying temp
    grid = _make_synthetic_grid(
        "single_monh_grid", teffs=[4800.0, 5000.0, 5200.0], loggs=[3.6, 4.0, 4.4],
        monhs=[0.0], temp_value=5000.0,
    )
    teff, logg = np.asarray(grid["teff"]), np.asarray(grid["logg"])
    grid["temp"] = _log_linear_temp(teff, logg)[:, None]
    return grid


def _log_linear_temp(teff, logg):
    # Linear in log10(temp), so the RBF's linear polynomial reproduces it exactly
    return 10.0 ** (3.7 + 1e-4 * (teff - 5000.0) + 0.02 * logg)


def test_rbf_handles_single_valued_metallicity_axis():
    # The single-valued [M/H] axis is left out of the fit: no crash, no warnings
    interpolator = AtmosphereInterpolator(interp="RBF")
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        atmo = interpolator.interp_atmo_grid(_single_monh_grid(), 4900.0, 3.8, 0.0)
    assert list(interpolator.rbf_grid.active_axes) == [True, True, False]
    assert np.all(np.isfinite(atmo.temp))
    assert np.allclose(atmo.temp, _log_linear_temp(4900.0, 3.8), rtol=1e-6)


def test_rbf_ignores_off_grid_value_on_single_valued_axis(caplog):
    # The grid says nothing about other [M/H], so the on-grid result is returned
    # and the request is reported as outside the grid
    grid = _single_monh_grid()
    interpolator = AtmosphereInterpolator(interp="RBF")
    on_grid = interpolator.interp_atmo_grid(grid, 4900.0, 3.8, 0.0)
    with caplog.at_level("INFO", logger="pysme.atmosphere.interpolation"):
        off_grid = interpolator.interp_atmo_grid(grid, 4900.0, 3.8, 0.3)

    assert np.array_equal(off_grid.temp, on_grid.temp)
    messages = [r.getMessage() for r in caplog.records]
    assert any("[M/H]=0.300 outside grid range" in m for m in messages)
    with pytest.raises(AtmosphereError, match="outside grid range"):
        interpolator.interp_atmo_grid(
            grid, 4900.0, 3.8, 0.3, interpolation_policy="error"
        )


def test_rbf_preserves_opflag_and_wlstd():
    # opflag and wlstd differ from the Atmo() defaults ([1]*20, 5000.0)
    teffs, loggs, monhs = [4900.0, 5100.0], [3.8, 4.2], [-0.2, 0.2]
    grid = _make_synthetic_grid(
        "metadata_grid", teffs, loggs, monhs, temp_value=5000.0,
        opflag_value=0, wlstd_value=6000.0,
    )
    interpolator = AtmosphereInterpolator(interp="RBF")
    atmo = interpolator.interp_atmo_grid(grid, 5000.0, 4.0, 0.0)

    assert np.array_equal(atmo.opflag, np.zeros(20, dtype=atmo.opflag.dtype))
    assert atmo.wlstd == 6000.0


def test_rbf_takes_opflag_and_wlstd_from_the_nearest_model_regardless_of_row_order():
    # RbfGrid sorts its grid points, so with reversed rows and a different
    # wlstd/opflag per model, a mix-up between the two orders shows
    teffs, loggs, monhs = [4900.0, 5100.0], [3.8, 4.2], [-0.2, 0.2]
    grid = _make_synthetic_grid(
        "row_order_grid", teffs, loggs, monhs, temp_value=5000.0
    )
    grid["wlstd"] = 1000.0 + np.arange(len(grid))
    grid["opflag"][:, 0] = np.arange(len(grid))
    grid = grid[np.arange(len(grid))[::-1]]

    query = (5100.0, 3.8, 0.2)  # on a grid point, so that model is the nearest
    expected = grid.get(*query)
    interpolator = AtmosphereInterpolator(interp="RBF")
    atmo = interpolator.interp_atmo_grid(grid, *query)

    assert atmo.wlstd == expected.wlstd
    assert np.array_equal(atmo.opflag, expected.opflag)


def test_rbf_reproduces_grid_models_on_any_scaled_value():
    # Teff steps of 100 and 250 K put Teff=4450 at 400.5 scaled units, and
    # logg=3.8 is stored as float32; the fit must still pass through these models
    teffs = [4000.0, 4100.0, 4200.0, 4450.0, 4700.0]
    loggs, monhs = [3.4, 3.8, 4.2], [0.0, 0.5]
    grid = _make_synthetic_grid(
        "uneven_teff_grid", teffs, loggs, monhs, temp_value=5000.0
    )
    rng = np.random.default_rng(0)
    grid["temp"] = rng.uniform(4000.0, 6000.0, size=grid["temp"].shape)

    query = (4450.0, 3.8, 0.5)
    expected = grid.get(*query)
    interpolator = AtmosphereInterpolator(interp="RBF")
    atmo = interpolator.interp_atmo_grid(grid, *query)

    assert np.allclose(atmo.temp, expected.temp, rtol=1e-9)


def test_rbf_index_space_maps_grid_steps_to_integers():
    axis_values = np.array([-5.0, -4.0, -3.0, -2.5, -2.0])
    values = [-6.0, -5.0, -4.5, -2.75, -2.0, -1.5]
    expected = [-1.0, 0.0, 0.5, 2.5, 4.0, 5.0]  # linear beyond the edges
    assert np.allclose(RbfGrid._to_index_space(values, axis_values), expected)
    assert np.array_equal(RbfGrid._to_index_space([0.3], np.array([0.0])), [0.0])


def _uneven_monh_grid():
    # Dense Teff and logg, with 1 dex [M/H] steps at the metal-poor end like MARCS
    teffs = np.arange(4000.0, 6001.0, 100.0)
    loggs = np.arange(1.0, 4.01, 0.5)
    monhs = [-5.0, -4.0, -3.0, -2.0, -1.5, -1.0, -0.75, -0.5, -0.25, 0.0]
    return _make_synthetic_grid(
        "uneven_monh_grid", teffs, loggs, monhs, temp_value=5000.0
    )


def test_rbf_index_selection_spreads_neighbors_over_uneven_metallicity_steps(
    monkeypatch,
):
    # In physical units the 1 dex steps are far apart, so the neighbors cover
    # only 2 [M/H] values; counting grid steps gives 4
    grid = _uneven_monh_grid()
    query = (5000.0, 2.5, -4.5)

    def n_monh(rbf_grid):
        neighbors = rbf_grid.neighbor_indices(query)
        return len(np.unique(rbf_grid.grid_coords[neighbors, 2]))

    index_grid = RbfGrid(grid, "PP")
    monkeypatch.setattr(RbfGrid, "_NEIGHBOR_SELECTION", "physical")
    physical_grid = RbfGrid(grid, "PP")

    assert n_monh(index_grid) == 4
    assert n_monh(physical_grid) == 2


@pytest.mark.parametrize("selection", ["index", "physical"])
def test_rbf_includes_all_neighbors_tied_with_the_last(monkeypatch, selection):
    # On a regular grid the neighbors of the central point are symmetric around
    # it; picking only some of the points tied at the cut would break that.
    # logg=3.8 is stored as float32, which must not break the ties either.
    monkeypatch.setattr(RbfGrid, "_NEIGHBOR_SELECTION", selection)
    steps = np.array([200.0, 0.2, 0.2])
    grid = _make_synthetic_grid(
        "regular_grid",
        teffs=np.arange(4000.0, 5601.0, 200.0),
        loggs=np.round(np.arange(3.0, 4.61, 0.2), 1),
        monhs=np.round(np.arange(-0.8, 0.81, 0.2), 1),
        temp_value=5000.0,
    )
    rbf_grid = RbfGrid(grid, "PP")
    query = np.array([4800.0, 3.8, 0.0])

    neighbors = rbf_grid.neighbor_indices(query)
    offsets = np.round((rbf_grid.grid_coords[neighbors] - query) / steps)
    offsets = {tuple(offset) for offset in offsets.astype(int)}

    assert len(neighbors) > RbfGrid._RBF_NEIGHBORS
    assert offsets == {tuple(-np.array(offset)) for offset in offsets}


def test_rbf_rejects_unknown_neighbor_selection(monkeypatch):
    monkeypatch.setattr(RbfGrid, "_NEIGHBOR_SELECTION", "nearest")
    with pytest.raises(AtmosphereError, match="unknown neighbor selection"):
        RbfGrid(_uneven_monh_grid(), "PP")


def test_rbf_drops_and_reports_model_with_non_positive_values(caplog):
    # Like the model in marcs2014.sav with a negative xna at its top depth point
    teffs, loggs, monhs = [4900.0, 5100.0], [3.8, 4.2], [-0.2, 0.2]
    grid = _make_synthetic_grid(
        "bad_model_grid", teffs, loggs, monhs, temp_value=5000.0
    )
    bad_row = 5  # (5100, 3.8, 0.2)
    grid["xna"][bad_row, 0] = -1.0

    interpolator = AtmosphereInterpolator(interp="RBF")
    with caplog.at_level("WARNING", logger="pysme.atmosphere.interpolation"):
        atmo = interpolator.interp_atmo_grid(grid, 5000.0, 4.0, 0.0)

    warnings = [r.getMessage() for r in caplog.records if r.levelname == "WARNING"]
    assert len(warnings) == 1
    assert "Teff=5100, log(g)=3.8, [M/H]=0.2" in warnings[0]
    assert "xna=-1 at depth point 0" in warnings[0]
    assert interpolator.rbf_grid.gridpoints == len(grid) - 1
    assert np.all(np.isfinite(atmo.xna))


def test_rbf_rejects_grid_with_several_models_at_one_point():
    teffs, loggs, monhs = [4900.0, 5100.0], [3.8, 4.2], [-0.2, 0.2]
    grid = _make_synthetic_grid(
        "duplicate_grid", teffs, loggs, monhs, temp_value=5000.0
    )
    grid["monh"][1] = grid["monh"][0]  # rows 0 and 1 now share (teff, logg, monh)

    interpolator = AtmosphereInterpolator(interp="RBF")
    with pytest.raises(AtmosphereError, match="one model per"):
        interpolator.interp_atmo_grid(grid, 5000.0, 4.0, 0.0)


def test_rbf_abundance_stays_consistent_with_interpolated_metallicity():
    # Each model's abundances are solar shifted by its own monh, like the
    # realized patterns in marcs2014.sav
    teffs, loggs, monhs = [4900.0, 5100.0], [3.8, 4.2], [-2.0, -1.0, 0.0]
    grid = _make_synthetic_grid(
        "abundance_grid", teffs, loggs, monhs, temp_value=5000.0,
        abund_shift=lambda monh: monh,
    )
    interpolator = AtmosphereInterpolator(interp="RBF")

    query_monh = -1.4  # between two grid values
    atmo = interpolator.interp_atmo_grid(grid, 5000.0, 4.0, query_monh)

    # The abundances must match the atmosphere's own interpolated monh (which is
    # only close to the query), not those of the nearest grid point
    expected_fe = _SOLAR_H12[elements_dict["Fe"]] + atmo.monh

    assert np.isclose(atmo.abund.get_pattern_abundance("Fe"), expected_fe, atol=1e-6)


def test_rbf_interpolates_spherical_geometry():
    # The only synthetic SPH test: height is interpolated together with radius
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
    # The helper uses radius=1e11 for SPH grids (1 for PP)
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

    # Default policy "allow": extrapolates instead of raising
    atmo = interpolator.interp_atmo_grid(grid, 50000.0, 4.0, 0.0)
    assert np.all(np.isfinite(atmo.temp))


def _real_depth_count(atmo):
    """
    Number of real depth points of a model. marcs2014.sav pads its depth arrays
    (to 72) with zeros or NaN, which the RBF output doesn't include.
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
    Spherical models of marcs2014.sav, a real grid with several metallicities,
    for both the structure and the abundances. Only checks that RBF reproduces
    a grid point and stays plausible between grid points, not that the
    abundance values themselves are physically correct.
    """
    grid = SavFile(lfs_atmo.get("marcs2014.sav"), source="marcs2014.sav", lfs=lfs_atmo)
    interpolator = AtmosphereInterpolator(interp="RBF", geom="SPH", lfs_atmo=lfs_atmo)
    logg, monh = 2.5, -1.0

    # At a grid point the RBF (no smoothing) reproduces the stored model
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

    # Between two grid points: finite, hotter with depth, and close to both
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
    Plane-parallel version of test_rbf_interpolation_matches_real_grid, on
    marcs2012p_t1.0.sav (56 depth points, no padding, 15 metallicities).
    """
    grid = SavFile(
        lfs_atmo.get("marcs2012p_t1.0.sav"), source="marcs2012p_t1.0.sav", lfs=lfs_atmo
    )
    interpolator = AtmosphereInterpolator(interp="RBF", geom="PP", lfs_atmo=lfs_atmo)
    logg, monh = 4.0, -1.0

    # At a grid point the RBF (no smoothing) reproduces the stored model
    teff = 5000.0
    atmo_interp = interpolator.interp_atmo_grid(grid, teff, logg, monh)
    atmo_grid = grid.get(teff, logg, monh)

    assert np.allclose(atmo_interp.temp, atmo_grid.temp, rtol=1e-6)
    assert np.isclose(
        atmo_interp.abund.get_pattern_abundance("Si"),
        atmo_grid.abund.get_pattern_abundance("Si"),
        rtol=1e-6,
    )

    # Between two grid points: finite, hotter with depth, and close to both
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
