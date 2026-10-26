"""Non-native regression tests for split fitting across automatic geometry."""
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from pysme.atmosphere.atmosphere import Atmosphere
from pysme.linelist.linelist import LineList
from pysme.sme import SME_Structure
from pysme.solve import SME_Solver, _least_squares_cost


def make_sme(*, geom=None, logg=3.2):
    sme = SME_Structure()
    sme.teff = 5250.0
    sme.logg = logg
    sme.monh = 0.0
    sme.atmo = Atmosphere(
        source="mock-geometry.sav",
        method="grid",
        geom=geom,
        depth="RHOX",
        interp="RHOX",
    )
    sme.wave = [np.array([5000.0, 5000.5, 5001.0])]
    sme.spec = [np.ones(3)]
    sme.uncs = [np.ones(3)]
    sme.mask = [np.ones(3, dtype=np.int16)]
    sme.linelist = LineList(
        pd.DataFrame(
            {
                "species": ["Fe 1"],
                "wlcent": [5000.5],
                "excit": [1.0],
                "gflog": [-1.0],
                "gamrad": [0.0],
                "gamqst": [0.0],
                "gamvw": [0.0],
                "atom_number": [26.0],
                "ionization": [1.0],
            }
        )
    )
    return sme


def make_solver(monkeypatch, outcomes, *, geometry_types_by_source=None):
    """Make an SME solver whose optimizer and synthesis are both local mocks."""
    import importlib

    solve_module = importlib.import_module("pysme.solve")
    if geometry_types_by_source is None:
        geometry_types_by_source = {
            "mock-geometry.sav": frozenset(("PP", "SPH")),
            "pure-pp.sav": frozenset(("PP",)),
            "pure-sph.sav": frozenset(("SPH",)),
        }
    calls = []
    synth_calls = []

    def fake_synthesize(local_sme, **_kwargs):
        requested = local_sme.atmo.geom_request
        source_geometries = geometry_types_by_source.get(
            local_sme.atmo.source, frozenset(("PP", "SPH"))
        )
        resolved = requested or (
            next(iter(source_geometries))
            if len(source_geometries) == 1
            else ("SPH" if local_sme.logg < 3.5 else "PP")
        )
        local_sme.atmo._set_resolved_geometry(resolved, requested)
        outcome = outcomes.get(resolved, outcomes.get("single"))
        flux = np.ones(3)
        if not isinstance(outcome, BaseException) and len(outcome) >= 3:
            flux = np.asarray(outcome[2], dtype=float)
        local_sme.synth = [flux]
        synth_calls.append(
            (id(local_sme), float(local_sme.logg), resolved, requested, flux.copy())
        )
        return local_sme

    solver = SME_Solver.__new__(SME_Solver)
    solver.synthesizer = SimpleNamespace(synthesize_spectrum=fake_synthesize)
    solver.filename = None
    solver.restore = False
    solver.parameter_names = []
    solver.derived_param = None
    solver.update_linelist = False
    solver._prepare_abundance_synthesis = False
    solver._smelib_state_prepared = False
    solver._latest_jacobian = np.zeros((3, 1))
    solver._can_prepare_abundance_synthesis = lambda *_args, **_kwargs: False
    solver._validate_field_against_wave = lambda *_args, **_kwargs: None
    solver._get_grid_geometries = lambda local_sme: geometry_types_by_source.get(
        local_sme.atmo.source, frozenset(("PP", "SPH"))
    )
    def fake_update_fitresults(local_sme, result, segments):
        if segments == "all":
            segments = list(range(len(local_sme.spec)))
        if local_sme.synth is None:
            local_sme.synth = [np.ones_like(np.asarray(values)) for values in local_sme.spec]
        residuals = []
        for segment in segments:
            mask = np.asarray(local_sme.mask_line[segment], dtype=bool)
            residuals.extend(
                ((np.asarray(local_sme.spec[segment]) - np.asarray(local_sme.synth[segment]))
                 / np.asarray(local_sme.uncs[segment]))[mask]
            )
            if local_sme.telluric is not None:
                residuals[-np.count_nonzero(mask):] = (
                    (np.asarray(local_sme.synth[segment])
                     * np.asarray(local_sme.telluric[segment])
                     - np.asarray(local_sme.spec[segment]))
                    / np.asarray(local_sme.uncs[segment])
                )[mask]
        local_sme.fitresults.residuals = np.asarray(residuals, dtype=float)
        local_sme.fitresults.chisq = float(
            np.sum(local_sme.fitresults.residuals**2)
        )
        local_sme.fitresults.values = np.asarray(result.x)
        local_sme.fitresults.parameters = list(solver.parameter_names)
        local_sme.fitresults.fit_uncertainties = np.zeros(len(solver.parameter_names))
        return local_sme

    solver.update_fitresults = fake_update_fitresults

    def fake_least_squares(*_args, **kwargs):
        bounds = np.asarray(kwargs["bounds"], dtype=float)
        residual_args = kwargs.get("args", ())
        local_sme = residual_args[0] if residual_args else None
        source_geometries = geometry_types_by_source.get(
            getattr(getattr(local_sme, "atmo", None), "source", None),
            frozenset(("PP", "SPH")),
        )
        if bounds.ndim == 2 and bounds.shape[0] == 2 and bounds.shape[1] == 1:
            lower, upper = bounds[:, 0]
            if len(source_geometries) == 1:
                side = next(iter(source_geometries))
            elif upper <= 3.5 + 1e-10:
                side = "SPH"
            elif lower >= 3.5 - 1e-10:
                side = "PP"
            else:
                side = ("SPH", "PP")[len(calls)]
        else:
            side = "single"
        calls.append((side, bounds.copy()))
        outcome = outcomes[side] if side in outcomes else outcomes["single"]
        if isinstance(outcome, BaseException):
            raise outcome
        value, cost = outcome[:2]
        return SimpleNamespace(
            x=np.array([value]),
            cost=float(cost),
            fun=np.array([np.sqrt(2 * float(cost))]),
            jac=np.ones((1, 1)),
            grad=np.zeros(1),
            optimality=0.0,
            active_mask=np.zeros(1, dtype=int),
            nfev=1,
            njev=1,
            status=1,
            message="mock success",
            success=True,
        )

    monkeypatch.setattr(solve_module, "least_squares", fake_least_squares)
    return solver, calls, synth_calls


def test_crossing_auto_fit_runs_both_sides_and_returns_best(monkeypatch):
    sme = make_sme()
    solver, calls, synth_calls = make_solver(
        monkeypatch,
        {
            "SPH": (3.4, 4.0, [1.5, 1.5, 1.5]),
            "PP": (3.6, 1.0, [1.0, 1.0, 1.0]),
            "single": (3.4, 4.0, [1.5, 1.5, 1.5]),
        },
    )

    result = solver.solve(
        sme,
        param_names=["logg"],
        bounds=np.array([[3.0], [4.0]]),
    )

    assert result is sme
    assert result.logg == pytest.approx(3.6)
    assert result.atmo.geom == "PP"
    assert result.atmo.geom_request is None
    assert {side for side, _ in calls} == {"SPH", "PP"}
    np.testing.assert_allclose(result.synth[0], [1.0, 1.0, 1.0])
    np.testing.assert_allclose(result.fitresults.residuals, [0.0, 0.0, 0.0])


@pytest.mark.parametrize(
    "source,expected_geom,accepted_logg",
    [
        ("pure-pp.sav", "PP", 3.2),
        ("pure-sph.sav", "SPH", 3.8),
    ],
)
def test_auto_geometry_pure_source_fits_once_across_boundary(
    monkeypatch, source, expected_geom, accepted_logg
):
    sme = make_sme()
    sme.atmo.source = source
    solver, calls, _ = make_solver(
        monkeypatch,
        {
            "SPH": (3.8, 1.0),
            "PP": (3.2, 1.0),
            "single": (accepted_logg, 1.0),
        },
    )

    result = solver.solve(
        sme,
        param_names=["logg"],
        bounds=np.array([[3.0], [4.0]]),
    )

    assert result is sme
    assert len(calls) == 1
    assert calls[0][0] == expected_geom
    assert result.logg == pytest.approx(accepted_logg)
    assert result.atmo.geom == expected_geom
    assert result.atmo.geom_request is None
    assert result.atmo.source == source


def test_crossing_auto_fit_exact_tie_prefers_sph(monkeypatch):
    sme = make_sme()
    solver, calls, _ = make_solver(
        monkeypatch,
        {"SPH": (3.4, 2.0), "PP": (3.6, 2.0), "single": (3.4, 2.0)},
    )

    result = solver.solve(
        sme,
        param_names=["logg"],
        bounds=np.array([[3.0], [4.0]]),
    )

    assert result is sme
    assert result.logg == pytest.approx(3.4)
    assert result.atmo.geom == "SPH"
    assert len(calls) == 2


@pytest.mark.parametrize("failed_side", ["SPH", "PP"])
def test_crossing_auto_fit_returns_successful_side_when_other_fails(
    monkeypatch, failed_side
):
    sme = make_sme()
    successful_side = "PP" if failed_side == "SPH" else "SPH"
    successful_value = 3.6 if successful_side == "PP" else 3.4
    solver, calls, _ = make_solver(
        monkeypatch,
        {
            "SPH": RuntimeError("mock SPH failure") if failed_side == "SPH" else (3.4, 2.0),
            "PP": RuntimeError("mock PP failure") if failed_side == "PP" else (3.6, 1.0),
            "single": (successful_value, 1.0),
        },
    )

    result = solver.solve(
        sme,
        param_names=["logg"],
        bounds=np.array([[3.0], [4.0]]),
    )

    assert result is sme
    assert result.logg == pytest.approx(successful_value)
    assert len(calls) == 2


def test_crossing_auto_fit_raises_when_both_sides_fail(monkeypatch):
    sme = make_sme()
    solver, calls, _ = make_solver(
        monkeypatch,
        {
            "SPH": RuntimeError("mock SPH failure"),
            "PP": RuntimeError("mock PP failure"),
            "single": RuntimeError("mock failure"),
        },
    )

    with pytest.raises(Exception):
        solver.solve(
            sme,
            param_names=["logg"],
            bounds=np.array([[3.0], [4.0]]),
        )
    assert len(calls) == 2


def test_crossing_auto_fit_uses_final_robust_objective(monkeypatch):
    """Final residuals can rank candidates differently from optimizer ``res.cost``."""
    sme = make_sme()
    sme.spec = [np.array([5.0, 5.0, 5.0])]
    sme.leastsquares_loss = "soft_l1"
    solver, calls, _ = make_solver(
        monkeypatch,
        {
            # SPH has lower mock optimizer cost, but higher soft-L1 cost.
            "SPH": (3.4, 1.0, [3.6, 3.6, 5.0]),
            "PP": (3.6, 2.0, [5.0, 3.0, 5.0]),
            "single": (3.4, 1.0, [3.6, 3.6, 5.0]),
        },
    )

    result = solver.solve(
        sme,
        param_names=["logg"],
        bounds=np.array([[3.0], [4.0]]),
    )

    assert result.logg == pytest.approx(3.6)
    assert {side for side, _ in calls} == {"SPH", "PP"}
    np.testing.assert_allclose(result.synth[0], [5.0, 3.0, 5.0])
    np.testing.assert_allclose(result.fitresults.residuals, [0.0, 2.0, 0.0])


def test_crossing_auto_fit_includes_telluric_in_final_objective(monkeypatch):
    sme = make_sme()
    sme.telluric = [np.array([0.5, 1.0, 2.0])]
    solver, calls, _ = make_solver(
        monkeypatch,
        {
            # Without telluric, PP has the smaller SSE. With telluric, SPH wins.
            "SPH": (3.4, 2.0, [0.8, 1.0, 1.0]),
            "PP": (3.6, 1.0, [1.0, 1.0, 1.1]),
            "single": (3.4, 2.0, [0.8, 1.0, 1.0]),
        },
    )

    result = solver.solve(
        sme,
        param_names=["logg"],
        bounds=np.array([[3.0], [4.0]]),
    )

    assert result.logg == pytest.approx(3.4)
    assert {side for side, _ in calls} == {"SPH", "PP"}


def make_region_solver(monkeypatch, outcomes, *, filename=None, restore=False):
    """Exercise _solve_geometry_regions with a mocked recursive solve call."""
    calls = []
    save_calls = []
    input_states = {}

    def fake_candidate_solve(candidate, local_sme, bounds, **_kwargs):
        gravity_column = 1 if np.asarray(bounds).shape[1] > 1 else 0
        lower, upper = np.asarray(bounds)[:, gravity_column]
        side = "SPH" if upper < 3.5 else "PP"
        calls.append((side, np.array(bounds, copy=True), candidate.filename, candidate.restore))
        input_states[side] = (local_sme.teff, local_sme.logg, local_sme.atmo.geom_request)
        outcome = outcomes[side]
        if isinstance(outcome, BaseException):
            local_sme.teff = 9999.0
            local_sme.atmo.geom = "SPH"
            raise outcome
        logg, cost = outcome
        local_sme.logg = logg
        local_sme.atmo._set_resolved_geometry("SPH" if side == "SPH" else "PP", None)
        local_sme.synth = [np.full(3, 0.8 if side == "SPH" else 1.0)]
        candidate._fit_cost = cost
        return local_sme

    monkeypatch.setattr(SME_Solver, "solve", fake_candidate_solve)
    monkeypatch.setattr(
        SME_Structure,
        "save",
        lambda local_sme, path: save_calls.append((local_sme, path)),
    )

    solver = SME_Solver.__new__(SME_Solver)
    solver.synthesizer = SimpleNamespace(
        wint={"stale": True}, known_sme=object(), atmosphere_interpolator=object()
    )
    solver.filename = filename
    solver.restore = restore
    solver.parameter_names = ["teff", "logg", "monh"]
    solver._fit_cost = None
    solver._latest_residual = None
    solver._latest_jacobian = None
    solver._smelib_state_prepared = False
    solver.progressbar = solver.progressbar_jacobian = None
    return solver, calls, save_calls, input_states


def test_geometry_regions_only_split_logg_and_preserve_input_bounds(monkeypatch):
    sme = make_sme()
    bounds = np.array([[5000.0, 3.0, -2.0], [6000.0, 4.0, 1.0]])
    original_bounds = bounds.copy()
    solver, calls, _, _ = make_region_solver(
        monkeypatch, {"SPH": (3.4, 2.0), "PP": (3.6, 1.0)}
    )

    result = solver._solve_geometry_regions(
        sme,
        bounds,
        param_names=solver.parameter_names,
        segments=[0],
    )

    assert result is sme
    np.testing.assert_array_equal(bounds, original_bounds)
    assert len(calls) == 2
    for side, region_bounds, _, _ in calls:
        np.testing.assert_array_equal(region_bounds[:, 0], original_bounds[:, 0])
        np.testing.assert_array_equal(region_bounds[:, 2], original_bounds[:, 2])
        if side == "SPH":
            assert region_bounds[0, 1] == original_bounds[0, 1]
            assert region_bounds[1, 1] < 3.5
        else:
            assert region_bounds[0, 1] == 3.5
            assert region_bounds[1, 1] == original_bounds[1, 1]


def test_geometry_regions_disables_candidate_persistence_and_saves_winner_once(
    monkeypatch, tmp_path
):
    sme = make_sme()
    filename = str(tmp_path / "winner.sme")
    solver, calls, save_calls, _ = make_region_solver(
        monkeypatch,
        {"SPH": (3.4, 2.0), "PP": (3.6, 1.0)},
        filename=filename,
        restore=True,
    )
    solver.parameter_names = ["logg"]

    result = solver._solve_geometry_regions(
        sme,
        np.array([[3.0], [4.0]]),
        param_names=["logg"],
        segments=[0],
    )

    assert result is sme
    assert all(candidate_filename is None for _, _, candidate_filename, _ in calls)
    assert all(candidate_restore is False for _, _, _, candidate_restore in calls)
    assert solver.filename == filename
    assert solver.restore is True
    assert len(save_calls) == 1
    assert save_calls[0][0] is result
    assert save_calls[0][1] == filename


def test_geometry_regions_skips_nonfinite_candidate_and_keeps_other(monkeypatch):
    sme = make_sme()
    solver, calls, _, _ = make_region_solver(
        monkeypatch, {"SPH": (3.4, np.nan), "PP": (3.6, 1.0)}
    )
    solver.parameter_names = ["logg"]

    result = solver._solve_geometry_regions(
        sme,
        np.array([[3.0], [4.0]]),
        param_names=["logg"],
        segments=[0],
    )

    assert result.logg == pytest.approx(3.6)
    assert len(calls) == 2


def test_geometry_regions_isolates_failed_candidate_state(monkeypatch):
    sme = make_sme()
    solver, calls, _, input_states = make_region_solver(
        monkeypatch,
        {"SPH": RuntimeError("SPH failed"), "PP": (3.6, 1.0)},
    )
    solver.parameter_names = ["logg"]

    result = solver._solve_geometry_regions(
        sme,
        np.array([[3.0], [4.0]]),
        param_names=["logg"],
        segments=[0],
    )

    assert result.logg == pytest.approx(3.6)
    assert input_states["SPH"][:2] == (5250.0, 3.2)
    assert input_states["PP"][:2] == (5250.0, 3.2)
    assert len(calls) == 2


@pytest.mark.parametrize("loss", ["linear", "soft_l1", "huber", "cauchy", "arctan"])
def test_least_squares_cost_matches_scipy_loss_definition(loss):
    residuals = np.array([0.5, -2.0])
    f_scale = 1.3
    z = (residuals / f_scale) ** 2
    if loss == "linear":
        rho = z
    elif loss == "soft_l1":
        rho = 2 * (np.sqrt(1 + z) - 1)
    elif loss == "huber":
        rho = np.where(z <= 1, z, 2 * np.sqrt(z) - 1)
    elif loss == "cauchy":
        rho = np.log1p(z)
    else:
        rho = np.arctan(z)
    expected = 0.5 * f_scale**2 * np.sum(rho)

    assert _least_squares_cost(residuals, loss, f_scale) == pytest.approx(expected)


@pytest.mark.parametrize(
    "geom,param_names,bounds",
    [
        ("PP", ["logg"], np.array([[3.0], [4.0]])),
        (None, ["logg"], np.array([[3.0], [3.4]])),
        (None, ["teff"], np.array([[5000.0], [5500.0]])),
    ],
)
def test_auto_geometry_keeps_single_side_paths(monkeypatch, geom, param_names, bounds):
    sme = make_sme(geom=geom)
    solver, calls, _ = make_solver(
        monkeypatch,
        {"SPH": (3.4, 1.0), "PP": (3.6, 1.0), "single": (3.4, 1.0)},
    )

    result = solver.solve(sme, param_names=param_names, bounds=bounds)

    assert result is sme
    assert len(calls) == 1


def test_update_fitresults_uses_telluric_for_residuals_and_chisq():
    """The persisted fit statistics must use the telluric-adjusted spectrum."""
    sme = make_sme()
    sme.spec = [np.array([1.0, 1.0, 1.0])]
    sme.synth = [np.array([0.8, 0.9, 1.1])]
    sme.uncs = [np.ones(3)]
    sme.telluric = [np.array([0.5, 2.0, 1.0])]

    solver = SME_Solver.__new__(SME_Solver)
    solver.parameter_names = ["logg"]
    solver.iteration = 7
    result = SimpleNamespace(
        x=np.array([3.2]),
        jac=np.ones((3, 1)),
        grad=np.zeros(1),
        fun=np.ones(3),
    )

    updated = solver.update_fitresults(sme, result, [0])

    expected_residuals = np.array([0.6, -0.8, -0.1])
    assert updated is sme
    np.testing.assert_allclose(sme.fitresults.residuals, expected_residuals)
    assert sme.fitresults.chisq == pytest.approx(
        np.sum(expected_residuals**2) / (expected_residuals.size - 1)
    )


@pytest.mark.parametrize(
    "bounds",
    [
        np.array([[3.0], [3.5]]),
        np.array([[3.5], [4.0]]),
    ],
)
def test_auto_geometry_boundary_at_3_5_stays_single_side(monkeypatch, bounds):
    """An exact 3.5 endpoint must not create a zero-width split candidate."""
    sme = make_sme()
    solver, calls, _ = make_solver(
        monkeypatch,
        {
            "SPH": (3.4, 1.0),
            "PP": (3.6, 1.0),
            "single": (3.4, 1.0),
        },
    )

    result = solver.solve(sme, param_names=["logg"], bounds=bounds)

    assert result is sme
    assert len(calls) == 1
    np.testing.assert_array_equal(calls[0][1], bounds)
