"""Geometry selection across successive model-atmosphere requests.

These tests deliberately exercise the Python-side policy at the MARCS mixed
grid boundary.  They do not initialize SMElib, and therefore remain useful as
regression tests for the atmosphere request/resolution state independently of
native synthesis.
"""
from copy import deepcopy
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
from flex.extensions.bindata import MultipleDataExtension

from pysme.abund import Abund
from pysme.atmosphere.atmosphere import Atmosphere, AtmosphereError, AtmosphereGrid
from pysme.atmosphere import interpolation
from pysme.atmosphere.interpolation import R_sun, logg_sun
from pysme.linelist.linelist import LineList
from pysme.sme import SME_Structure
from pysme.solve import SME_Solver
from pysme.synthesize import Synthesizer


def make_grid(geometries):
    """Build a small MARCS-like grid with separate PP and SPH families."""
    teff_values = (5000.0, 5500.0)
    grid = AtmosphereGrid(len(teff_values) * len(geometries), 6)
    grid.source = "mixed-geometry.sav"
    grid.abund_format = "H=12"
    abund = Abund.solar().get_pattern(type="H=12", raw=True)
    for i, (teff, (logg, geom)) in enumerate(
        (t, model) for t in teff_values for model in geometries
    ):
        grid.teff[i] = teff
        grid.logg[i] = logg
        grid.monh[i] = 0.0
        grid.vturb[i] = 1.0
        grid.lonh[i] = 1.5
        is_spherical = geom == "SPH"
        grid.radius[i] = (
            R_sun * 10 ** ((logg_sun - logg) / 2) if is_spherical else 1.0
        )
        grid.height[i] = np.linspace(3e7, -1e7, 6) if is_spherical else 0.0
        grid.wlstd[i] = 5000.0
        grid.opflag[i] = 1
        grid.abund[i] = abund
        grid.rhox[i] = np.geomspace(1e-4, 1e-1, 6)
        grid.tau[i] = np.geomspace(1e-5, 1e-2, 6)
        grid.temp[i] = teff * np.linspace(0.8, 1.2, 6)
        grid.rho[i] = np.geomspace(1e-10, 1e-7, 6)
        grid.xna[i] = np.geomspace(1e14, 1e17, 6)
        grid.xne[i] = np.geomspace(1e10, 1e13, 6)
    return grid


@pytest.fixture
def mixed_geometry_grid(monkeypatch):
    grid = make_grid(
        ((2.0, "SPH"), (3.0, "SPH"), (3.5, "SPH"),
         (4.0, "PP"), (4.5, "PP"))
    )
    monkeypatch.setattr(interpolation, "SavFile", lambda *args, **kwargs: grid)
    return grid


@pytest.fixture
def pure_pp_grid(monkeypatch):
    grid = make_grid(
        ((2.0, "PP"), (3.0, "PP"), (3.5, "PP"),
         (4.0, "PP"), (4.5, "PP"))
    )
    grid.source = "pure-pp.sav"
    monkeypatch.setattr(interpolation, "SavFile", lambda *args, **kwargs: grid)
    return grid


def make_synthesizer():
    # Exercise the real atmosphere pipeline without initializing a native DLL.
    synth = Synthesizer.__new__(Synthesizer)
    synth.atmosphere_interpolator = None
    synth.lfs_atmo = SimpleNamespace(get=lambda name: name)
    return synth


def make_sme(geom=None):
    sme = SME_Structure()
    sme.teff = 5250.0
    sme.monh = 0.0
    sme.atmo = Atmosphere(
        source="mixed-geometry.sav", method="grid", geom=geom,
        depth="RHOX", interp="RHOX",
    )
    return sme


def test_empty_geometry_constructor_requests_auto():
    atmo = Atmosphere(geom="")
    assert atmo.geom_request is None


@pytest.mark.filterwarnings("ignore:Covariance of the parameters could not be estimated")
@pytest.mark.parametrize("gravities", [(2.5, 4.0, 2.5), (4.0, 2.5, 4.0)])
def test_auto_geometry_matches_fresh_requests(mixed_geometry_grid, gravities):
    synth = make_synthesizer()
    sme = make_sme()
    for logg in gravities:
        sme.logg = logg
        synth.get_atmosphere(sme)

        fresh = make_sme()
        fresh.logg = logg
        make_synthesizer().get_atmosphere(fresh)

        assert sme.atmo.geom == fresh.atmo.geom
        for field in ("temp", "rhox", "rho", "height"):
            np.testing.assert_allclose(sme.atmo[field], fresh.atmo[field])


@pytest.mark.filterwarnings("ignore:Covariance of the parameters could not be estimated")
@pytest.mark.parametrize(
    "logg, expected",
    [(3.499999, "SPH"), (3.5, "PP"), (3.500001, "PP")],
)
def test_auto_geometry_boundary_policy_and_request_state(
    mixed_geometry_grid, logg, expected
):
    sme = make_sme()
    sme.logg = logg
    make_synthesizer().get_atmosphere(sme)

    assert sme.atmo.geom_request is None
    assert sme.atmo.geom == expected


def test_auto_spherical_geometry_rejects_tau_transfer(mixed_geometry_grid):
    sme = make_sme()
    sme.atmo.depth = "TAU"
    sme.atmo.interp = "TAU"
    sme.logg = 2.5

    with pytest.raises(AtmosphereError, match="requires.*RHOX.*not supported"):
        make_synthesizer().get_atmosphere(sme)


def test_solver_final_synthesis_uses_accepted_auto_geometry(
    monkeypatch, mixed_geometry_grid
):
    """The returned model must match the accepted parameters, not the last trial."""
    import importlib

    solve_module = importlib.import_module("pysme.solve")
    sme = make_sme()
    sme.logg = 2.5
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

    captured = {}

    def fake_synthesize(local_sme, **_kwargs):
        captured["accepted_logg"] = local_sme.logg
        make_synthesizer().get_atmosphere(local_sme)
        captured["resolved_geom"] = local_sme.atmo.geom
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
    solver.get_bounds = lambda _sme: (np.array([2.0]), np.array([4.5]))
    solver.get_step_sizes = lambda _names: np.array([0.1])
    solver.get_default_values = lambda _sme: np.array([2.5])
    solver.update_fitresults = lambda local_sme, _res, _segments: local_sme

    def fake_least_squares(*_args, **_kwargs):
        return SimpleNamespace(x=np.array([3.75]), jac=np.zeros((3, 1)))

    monkeypatch.setattr(solve_module, "least_squares", fake_least_squares)
    result = solver.solve(sme, param_names=["logg"])

    assert result is sme
    assert captured["accepted_logg"] == pytest.approx(3.75)
    assert captured["resolved_geom"] == "PP"
    assert result.atmo.geom_request is None


@pytest.mark.filterwarnings("ignore:Covariance of the parameters could not be estimated")
def test_auto_geometry_survives_multiple_boundary_crossings(mixed_geometry_grid):
    synth = make_synthesizer()
    sme = make_sme()
    sequence = (2.5, 3.499999, 3.5, 3.500001, 4.25, 3.0, 2.5)

    for logg in sequence:
        sme.logg = logg
        synth.get_atmosphere(sme)

        fresh = make_sme()
        fresh.logg = logg
        make_synthesizer().get_atmosphere(fresh)

        assert sme.atmo.geom_request is None
        assert sme.atmo.geom == fresh.atmo.geom
        for field in ("temp", "rhox", "rho", "height"):
            np.testing.assert_allclose(sme.atmo[field], fresh.atmo[field])


@pytest.mark.filterwarnings("ignore:Covariance of the parameters could not be estimated")
@pytest.mark.parametrize("geom,gravities", [("PP", (4.0, 4.5, 4.25)), ("SPH", (2.5, 3.0, 3.25))])
def test_explicit_geometry_stays_fixed(mixed_geometry_grid, geom, gravities):
    synth = make_synthesizer()
    sme = make_sme(geom)
    for logg in gravities:
        sme.logg = logg
        synth.get_atmosphere(sme)
        assert sme.atmo.geom == geom
        assert sme.atmo.geom_request == geom


@pytest.mark.filterwarnings("ignore:Covariance of the parameters could not be estimated")
def test_explicit_override_can_be_returned_to_auto(mixed_geometry_grid):
    synth = make_synthesizer()
    sme = make_sme()

    sme.atmo.geom = "PP"
    sme.logg = 4.25
    synth.get_atmosphere(sme)
    assert sme.atmo.geom_request == "PP"
    assert sme.atmo.geom == "PP"

    sme.atmo.geom = None
    sme.logg = 2.5
    synth.get_atmosphere(sme)
    assert sme.atmo.geom_request is None
    assert sme.atmo.geom == "SPH"


@pytest.mark.filterwarnings("ignore:Covariance of the parameters could not be estimated")
def test_auto_state_survives_deepcopy_and_new_synthesizer(mixed_geometry_grid):
    sme = make_sme()
    sme.logg = 2.5
    make_synthesizer().get_atmosphere(sme)

    copied = deepcopy(sme)
    assert copied.atmo.geom_request is None
    copied.logg = 4.25
    make_synthesizer().get_atmosphere(copied)
    assert copied.atmo.geom_request is None
    assert copied.atmo.geom == "PP"


@pytest.mark.filterwarnings("ignore:Covariance of the parameters could not be estimated")
def test_auto_state_roundtrips_through_dict(mixed_geometry_grid):
    sme = make_sme()
    sme.logg = 2.5
    make_synthesizer().get_atmosphere(sme)
    assert sme.atmo.geom == "SPH"
    assert sme.atmo.geom_request is None

    restored_atmo = Atmosphere.from_dict(sme.atmo.to_dict())
    assert restored_atmo.geom == "SPH"
    assert restored_atmo.geom_request is None


@pytest.mark.filterwarnings("ignore:Covariance of the parameters could not be estimated")
def test_auto_state_roundtrips_through_sme(mixed_geometry_grid, tmp_path):
    sme = make_sme()
    sme.logg = 2.5
    make_synthesizer().get_atmosphere(sme)
    assert sme.atmo.geom == "SPH"
    assert sme.atmo.geom_request is None

    path = tmp_path / "auto-geometry.sme"
    sme.save(str(path))
    restored_sme = SME_Structure.load(str(path))
    assert restored_sme.atmo.geom == "SPH"
    assert restored_sme.atmo.geom_request is None


def test_legacy_atmosphere_header_treats_resolved_geometry_as_explicit(
    mixed_geometry_grid,
):
    sme = make_sme()
    sme.logg = 2.5
    make_synthesizer().get_atmosphere(sme)
    extension = sme.atmo._save()
    header = dict(extension.header)
    header.pop("geom_request", None)
    legacy = MultipleDataExtension(header, dict(extension.data))

    restored = Atmosphere._load(legacy)
    assert restored.geom == "SPH"
    assert restored.geom_request == "SPH"


@pytest.mark.filterwarnings("ignore:Covariance of the parameters could not be estimated")
def test_auto_pure_pp_grid_keeps_pp_below_mixed_boundary(pure_pp_grid):
    sme = make_sme()
    sme.atmo.source = "pure-pp.sav"
    sme.logg = 2.5
    make_synthesizer().get_atmosphere(sme)

    assert sme.atmo.geom_request is None
    assert sme.atmo.geom == "PP"
