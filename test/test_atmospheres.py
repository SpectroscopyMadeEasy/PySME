# -*- coding: utf-8 -*-
# TODO implement atmosphere tests
import numpy as np
import pytest

from pysme.atmosphere.atmosphere import Atmosphere
from pysme.atmosphere.interpolation import AtmosphereInterpolator, R_sun, logg_sun
from pysme.atmosphere.savfile import SavFile

from .test_largefilestorage import lfs_atmo, skipif_lfs

pytestmark = pytest.mark.filterwarnings(
    "ignore:Covariance of the parameters could not be estimated:scipy.optimize.OptimizeWarning"
)


@pytest.fixture
def atmosphere_name():
    # TODO iterate over all possible options
    return "marcs2012p_t1.0.sav"


@pytest.fixture
def atmosphere(atmosphere_name):
    atmo = Atmosphere(source=atmosphere_name, method="grid", interp="TAU", depth="RHOX")
    return atmo


@pytest.fixture
def interpolator(atmosphere, lfs_atmo):
    interp = AtmosphereInterpolator(
        depth=atmosphere.depth,
        interp=atmosphere.depth,
        geom=atmosphere.geom,
        lfs_atmo=lfs_atmo,
    )
    return interp


@pytest.fixture
def atmosphere_grid(atmosphere_name, lfs_atmo):
    name = lfs_atmo.get(atmosphere_name)
    atmo = SavFile(name, source=name)
    return atmo


@skipif_lfs
@pytest.mark.usefixtures("lfs_atmo")
def test_grid_point(atmosphere_name, atmosphere_grid, lfs_atmo, interpolator):
    # TODO: get this values from the grid
    teff = 7000
    logg = 4
    monh = 0

    atmo_interp = interpolator.interp_atmo_grid(atmosphere_name, teff, logg, monh)
    atmo_grid = atmosphere_grid.get(teff, logg, monh)

    assert np.allclose(atmo_interp.temp, atmo_grid.temp[1:])
    assert np.allclose(atmo_interp.tau, atmo_grid.tau[1:])
    assert np.allclose(atmo_interp.rhox, atmo_grid.rhox[1:])
    assert np.allclose(atmo_interp.rho, atmo_grid.rho[1:])
    assert np.allclose(atmo_interp.xna, atmo_grid.xna[1:])
    assert np.allclose(atmo_interp.xne, atmo_grid.xne[1:])


def _make_spherical_test_atmo(temp_offset=0.0, height_shift=0.0):
    ndep = 5
    atmo = Atmosphere(interp="RHOX")
    atmo.teff = 4500.0 + temp_offset
    atmo.logg = 2.0
    atmo.monh = -0.5
    atmo.vturb = 1.5
    atmo.lonh = 1.5
    # Radius consistent with logg=2.0 (a giant) at ~1 solar mass, via the same
    # mass-conservation relation the production code uses (see interp_atmo_pair).
    atmo.radius = R_sun * 10 ** ((logg_sun - atmo.logg) * 0.5)
    atmo.rhox = np.array([1e-4, 3e-4, 1e-3, 3e-3, 1e-2], dtype=float)
    atmo.tau = np.array([1e-5, 3e-5, 1e-4, 3e-4, 1e-3], dtype=float)
    atmo.temp = np.array([4000.0, 4200.0, 4400.0, 4600.0, 4800.0], dtype=float) + temp_offset
    atmo.xne = np.array([1e10, 2e10, 4e10, 8e10, 1.6e11], dtype=float)
    atmo.xna = np.array([1e14, 2e14, 4e14, 8e14, 1.6e15], dtype=float)
    atmo.rho = np.array([1e-10, 2e-10, 4e-10, 8e-10, 1.6e-9], dtype=float)
    atmo.height = np.array([-5e7, -3e7, -1e7, 1e7, 3e7], dtype=float) + height_shift
    return atmo


def test_interp_atmo_pair_interpolates_spherical_height():
    interpolator = AtmosphereInterpolator(depth="RHOX", interp="RHOX", geom="SPH")
    atmo1 = _make_spherical_test_atmo(temp_offset=0.0, height_shift=0.0)
    atmo2 = _make_spherical_test_atmo(temp_offset=200.0, height_shift=2e7)

    out = interpolator.interp_atmo_pair(atmo1, atmo2, frac=0.5, interpvar="RHOX")

    assert np.all(np.isfinite(out.height))
    assert len(out.height) == len(out.temp)
    assert not np.allclose(out.height, atmo1.height[: len(out.height)])
    assert not np.allclose(out.height, atmo2.height[: len(out.height)])

    # Coarse check: combining log(height+radius) is close to averaging height directly
    # when radius dominates. atol accounts for curve_fit's registration precision on
    # log10(height+radius)~12, which amplifies to an absolute noise floor of order
    # radius * (fit tolerance) ~ 1e4 cm once converted back out of log-space.
    assert np.allclose(
        out.height, 0.5 * (atmo1.height + atmo2.height), atol=1e5
    )
    # Exact check: at frac=0.5 with equal radius/logg, the combined-log result should be
    # the geometric mean of (height + radius), not the arithmetic mean used above; atol
    # allows for the same curve_fit registration noise described above.
    radius = atmo1.radius
    expected_height = np.sqrt((atmo1.height + radius) * (atmo2.height + radius)) - radius
    assert np.allclose(
        out.height, expected_height[: len(out.height)], atol=5e4, rtol=1e-3
    )
