"""Physical acceptance tests for the bundled MARCS 2012/2014 grid interfaces.

Small, real grid samples always run offline. CI additionally preloads the full
grids and sets PYSME_TEST_OFFICIAL_MARCS=1: load failures must then fail, not skip.
These checks apply to these standard-composition grids, not arbitrary user
abundances or custom MOD files. They inspect atmo.abund, never a replacement
sme.abund that could conceal an atmosphere-reader error.
"""

from copy import copy
import hashlib
import json
import os
from pathlib import Path

import numpy as np
import pytest

from pysme.atmosphere.interpolation import AtmosphereInterpolator
from pysme.atmosphere.savfile import SavFile
from pysme.large_file_storage import setup_atmo
from pysme.sme import SME_Structure


DATA = Path(__file__).with_name("data") / "official_marcs"
PROVENANCE = json.loads((DATA / "provenance.json").read_text())
GRID_NAMES = ("marcs2012.sav", "marcs2014.sav")
REFERENCE_POINTS = ((5750, 4.5, 0), (4500, 1.5, -2), (6000, 4, -2))
INTERPOLATION_POINTS = ((5500, 4, -2), (5625, 4.25, -1.75))


def assert_standard_h_he(atmosphere, context, expected_he=None):
    """Independent physical guard: these grids have He/H about 0.085."""
    h12 = atmosphere.abund(type="H=12")
    ratio = atmosphere.abund(type="n/nH")["He"]
    message = (
        f"{context}: Teff={atmosphere.teff}, logg={atmosphere.logg}, "
        f"[M/H]={atmosphere.monh}, A(H)={h12['H']}, "
        f"A(He)={h12['He']}, N_He/N_H={ratio}"
    )
    assert np.isfinite(h12["H"]) and np.isfinite(h12["He"]), message
    assert h12["H"] == pytest.approx(12.0, abs=1e-8), message
    assert np.isfinite(ratio) and 0.06 < ratio < 0.12, message
    if expected_he is not None:
        assert h12["He"] == pytest.approx(expected_he, abs=2e-6), message


@pytest.fixture(scope="module", params=GRID_NAMES)
def sample_grid(request):
    name = request.param
    metadata = PROVENANCE["grids"][name]
    path = DATA / metadata["fixture"]
    assert hashlib.sha256(path.read_bytes()).hexdigest() == metadata["fixture_sha256"]
    # lfs=None keeps checked-in fixtures read-only.
    return SavFile(path, source=name), metadata


@pytest.mark.parametrize("point", REFERENCE_POINTS)
def test_real_sample_grid_record_h_he(sample_grid, point):
    grid, metadata = sample_grid
    record = next(
        r for r in metadata["records"]
        if (r["teff"], r["logg"], r["monh"]) == point
    )
    atmosphere = grid.get(*point)
    # Golden numbers come from original linear He/H, not Abund.totype.
    assert_standard_h_he(atmosphere, grid.source, record["expected_A_He"])


@pytest.mark.parametrize("point", INTERPOLATION_POINTS)
@pytest.mark.filterwarnings("ignore:Covariance of the parameters could not be estimated")
def test_real_sample_interpolation_h_he(sample_grid, point):
    grid, _ = sample_grid
    interpolator = AtmosphereInterpolator(
        depth="RHOX", interp="RHOX", lfs_atmo=object()
    )
    atmosphere = interpolator.interp_atmo_grid(grid, *point)
    assert_standard_h_he(atmosphere, f"{grid.source}: interpolation", 10.927363695)


def assert_sme_roundtrip(atmosphere, path, context):
    expected = atmosphere.abund(type="H=12")["He"]
    cloned = copy(atmosphere.abund)
    assert cloned(type="H=12")["He"] == pytest.approx(expected, abs=2e-6)
    sme = SME_Structure()
    sme.atmo = atmosphere
    # Deliberately leave sme.abund alone; it is not the quantity under test.
    sme.save(str(path))
    restored = SME_Structure.load(str(path))
    assert_standard_h_he(restored.atmo, f"{context}: SME round trip", expected)


def test_real_sample_sme_persistence_h_he(sample_grid, tmp_path):
    grid, _ = sample_grid
    atmosphere = grid.get(4500, 1.5, -2)
    assert_standard_h_he(atmosphere, grid.source, 10.927363695)
    assert_sme_roundtrip(atmosphere, tmp_path / "sample.sme", grid.source)


@pytest.fixture(scope="module", params=GRID_NAMES)
def full_grid(request):
    if os.environ.get("PYSME_TEST_OFFICIAL_MARCS") != "1":
        pytest.skip("Full-grid acceptance is opt-in locally and mandatory in CI")
    name = request.param
    storage = setup_atmo()
    # No availability skip, exception catch or synthetic fallback here.
    path = storage.get(name)
    return SavFile(path, source=name)


def test_full_official_grid_h_he(full_grid):
    grid = full_grid
    assert len(grid) > 4000, f"{grid.source}: expected the full official grid"
    # Exercise actual Atmosphere construction for every grid record, including
    # all metallicities, rather than only decoding the raw abundance arrays.
    for index in range(len(grid)):
        atmosphere = grid[index]
        raw_h, raw_he = grid["abund"][index, :2]
        context = f"{grid.source}: record {index}, raw H/He={raw_h}/{raw_he}"
        assert_standard_h_he(atmosphere, context)


@pytest.mark.filterwarnings("ignore:Covariance of the parameters could not be estimated")
def test_full_official_grid_interpolation_and_persistence_h_he(full_grid, tmp_path):
    grid = full_grid
    for point in REFERENCE_POINTS:
        atmosphere = grid.get(*point)
        expected = 10.92783498 if point[2] == 0 else 10.927363695
        assert_standard_h_he(atmosphere, grid.source, expected)
    interpolator = AtmosphereInterpolator(
        depth="RHOX", interp="RHOX", lfs_atmo=object()
    )
    for index, point in enumerate(INTERPOLATION_POINTS):
        atmosphere = interpolator.interp_atmo_grid(grid, *point)
        assert_standard_h_he(atmosphere, f"{grid.source}: interpolation", 10.927363695)
        assert_sme_roundtrip(atmosphere, tmp_path / f"interpolated-{index}.sme", grid.source)


def test_standard_grid_guard_rejects_order_of_magnitude_helium_error(sample_grid):
    grid, _ = sample_grid
    atmosphere = grid.get(4500, 1.5, -2)
    # Ensure the acceptance criterion itself catches the historical-sized error.
    atmosphere.abund.set_A("He", 12.113)
    with pytest.raises(AssertionError, match="N_He/N_H"):
        assert_standard_h_he(atmosphere, "deliberately corrupted sample")
