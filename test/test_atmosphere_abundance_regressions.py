"""Regression tests for atmosphere abundance formats, persistence and readers."""

from pathlib import Path
from copy import copy

import numpy as np
import pytest

from pysme.abund import Abund
from pysme.atmosphere.atmosphere import Atmosphere, AtmosphereGrid
from pysme.atmosphere.krzfile import KrzFile, atmoic_mass
from pysme.atmosphere.savfile import SavFile


def legacy_pattern():
    pattern = np.full(99, -99.0)
    pattern[:2] = [0.922, 0.078]
    pattern[25] = -4.54
    return pattern


@pytest.mark.parametrize("monh", [-3.0, 0.0, 0.5])
def test_legacy_helium_and_metallicity(monh):
    abund = Abund(monh=monh, pattern=legacy_pattern(), type="sme_legacy")()
    assert abund["H"] == 12
    assert abund["He"] == pytest.approx(12 + np.log10(0.078 / 0.922))
    assert abund["Fe"] == pytest.approx(12 - 4.54 - np.log10(0.922) + monh)


@pytest.mark.parametrize("fmt", ["sme", "sme_legacy", "kurucz"])
@pytest.mark.parametrize("he", [10.0, 10.93, 12.0])
def test_abundance_roundtrip_and_independent_normalization(fmt, he):
    original = np.full(99, np.nan)
    original[[0, 1, 25]] = [12, he, 7.5]
    converted = Abund.totype(original, fmt, raw=True)
    he_over_h = 10 ** (he - 12)
    h_fraction = 1 / (1 + he_over_h + 10 ** (7.5 - 12))
    assert converted[0] == pytest.approx(h_fraction)
    if fmt == "kurucz":
        # PySME's documented Kurucz convention: metals relative to H+He.
        assert converted[25] == pytest.approx(-4.5 - np.log10(1 + he_over_h))
    else:
        assert converted[25] == pytest.approx(-4.5 + np.log10(h_fraction))
    expected_he = h_fraction * he_over_h
    if fmt == "sme":
        expected_he = np.log10(expected_he)
    assert converted[1] == pytest.approx(expected_he)
    np.testing.assert_allclose(Abund.fromtype(converted, fmt, raw=True), original)


def make_grid(fmt="sme"):
    grid = AtmosphereGrid(2, 3)
    for name in grid.dtype.names:
        grid[name] = 0
    grid["teff"] = [5000, 5500]
    grid["logg"] = 4
    grid["monh"] = [-2, 0]
    grid["temp"] = [4000, 4500, 5000]
    grid["abund"] = legacy_pattern()
    grid.abund_format = fmt
    return grid


@pytest.mark.parametrize("fmt", ["sme", "sme_legacy", "kurucz", "H=12"])
def test_grid_record_honors_format(fmt):
    pattern = Abund(pattern="grevesse2007", monh=0)(raw=True)
    grid = make_grid(fmt)
    grid["abund"] = Abund.totype(pattern, fmt, raw=True)
    expected = pattern.copy()
    expected[2:] -= 2
    np.testing.assert_allclose(grid[0].abund(raw=True), expected, atol=2e-6)
    np.testing.assert_allclose(grid[:1].abund(raw=True), expected, atol=2e-6)
    assert grid.ndep == 3


@pytest.fixture
def empty_sav_cache(monkeypatch):
    monkeypatch.setattr(SavFile, "_cache", {})


def test_old_numpy_cache_metadata_and_multiple_files(tmp_path, empty_sav_cache):
    first_path = tmp_path / "first.sav.npz"
    second_path = tmp_path / "second.sav.npz"
    legacy = make_grid()
    legacy.save(first_path)
    modern = make_grid()
    modern["abund"][:, 1] = np.log10(0.078)
    modern.save(second_path)
    first = SavFile(first_path)
    second = SavFile(second_path)
    assert first.abund_format == "sme_legacy"
    assert second.abund_format == "sme"
    expected = 12 + np.log10(0.078 / 0.922)
    assert first[0].abund()["He"] == pytest.approx(expected, abs=1e-6)
    assert second[0].abund()["He"] == pytest.approx(expected, abs=1e-6)
    np.testing.assert_array_equal(first["temp"], legacy["temp"])
    np.testing.assert_array_equal(first["abund"], legacy["abund"])
    assert SavFile(first_path) is first
    assert SavFile(str(second_path)) is second
    assert len(SavFile._cache) == 2
    # New metadata persists and a later load does not reinterpret it again.
    migrated_path = tmp_path / "migrated.npz"
    first.save(migrated_path)
    assert SavFile(migrated_path).abund_format == "sme_legacy"


def test_raw_idl_grid_legacy_metadata(monkeypatch, tmp_path, empty_sav_cache):
    grid = make_grid()
    raw = dict(atmo_grid_maxdep=3, atmo_grid_natmo=2,
               atmo_grid=grid.view(np.recarray))
    calls = []
    def fake_readsav(filename):
        calls.append(filename)
        return raw
    monkeypatch.setattr("pysme.atmosphere.savfile.readsav", fake_readsav)
    path = tmp_path / "marcs.sav"
    path.write_bytes(b"synthetic IDL input")
    result = SavFile(path)
    assert result.abund_format == "sme_legacy"
    assert result[0].abund()["He"] == pytest.approx(10.927363695, abs=1e-6)
    assert SavFile(path) is result
    assert len(calls) == 1


def test_mixed_helium_formats_rejected(tmp_path, empty_sav_cache):
    grid = make_grid()
    grid["abund"][1, 1] = -1.1
    path = tmp_path / "mixed.npz"
    grid.save(path)
    with pytest.raises(ValueError, match="inconsistent helium"):
        SavFile(path)


def test_marcs_krz_keeps_log_helium(tmp_path):
    pattern = legacy_pattern()
    pattern[1] = np.log10(pattern[1])
    values = np.append(pattern, 2).reshape(10, 10)
    lines = ["MARCS TEFF=5000 GRAVITY=4.0\n",
             "MODEL TYPE=0 WLSTD=5000 VTURB=2 L/H=1.25\n",
             " ".join(["1"] * 20) + " - opacity\n"]
    lines += [" ".join(map(str, row)) + "\n" for row in values]
    lines += ["0.1,4500,1e10,1e15,1e-9\n", "1.0,5000,2e10,2e15,2e-9\n"]
    path = tmp_path / "marcs.krz"
    path.write_text("".join(lines))
    result = KrzFile(path)
    assert result.abund()["He"] == pytest.approx(12 + np.log10(0.078 / 0.922))
    assert result.abund()["Fe"] == pytest.approx(12 - 4.54 - np.log10(0.922))


@pytest.mark.parametrize("monh", [-2.0, 0.0, 0.5])
def test_atlas_density_uses_atomic_pressure_and_effective_abundances(tmp_path, monh):
    text = Path(__file__).with_name("testatmo1.krz").read_text()
    text = text.replace("[0.0]", f"[{monh}]").replace(
        "ABUNDANCE SCALE   1.00000", f"ABUNDANCE SCALE   {10**monh:.8E}")
    path = tmp_path / "atlas.krz"
    path.write_text(text)
    result = KrzFile(path)
    expected_xna = result.P_gas / (1.380649e-16 * result.temp) - result.xne
    np.testing.assert_allclose(result.xna, expected_xna)
    ratios = result.abund(type="n/nH")
    total_mass = sum(ratios[e] * atmoic_mass[e] for e in ratios if np.isfinite(ratios[e]))
    total_number = sum(v for v in ratios.values() if np.isfinite(v))
    expected_rho = expected_xna * (total_mass / total_number) * 1.66053906660e-24
    np.testing.assert_allclose(result.rho, expected_rho)
    # Composition containing just H, He and Fe makes the expected mu explicit.
    pattern = np.full(99, np.nan)
    pattern[[0, 1, 25]] = [12, 10.93, 7.5]
    result.abund = Abund(monh=monh, pattern=pattern, type="H=12")
    he_ratio = 10**(10.93 - 12)
    fe_ratio = 10**(7.5 + monh - 12)
    expected_mu = (1.008 + 4.002602 * he_ratio + 55.845 * fe_ratio) / (1 + he_ratio + fe_ratio)
    assert result.get_mu_from_abund() == pytest.approx(expected_mu)


def test_atlas_header_decimals_signs_and_scientific_abundances(tmp_path):
    text = Path(__file__).with_name("testatmo1.krz").read_text()
    text = text.replace("GRAVITY 4.50000", "GRAVITY -0.50000")
    text = text.replace("VTURB=2", "VTURB=2.75")
    text = text.replace("1 0.92040 2 0.07834", "1 9.2040E-1 2 7.834E-2")
    path = tmp_path / "atlas.krz"
    path.write_text(text)
    result = KrzFile(path)
    assert result.logg == -0.5
    assert result.vturb == 2.75
    assert result.abund()["He"] == pytest.approx(12 + np.log10(0.07834 / 0.9204), abs=1e-6)


@pytest.mark.parametrize("fmt", ["sme", "sme_legacy", "kurucz", "H=12", "n/nTot", "n/nH"])
def test_flex_abundance_and_atmosphere_preserve_internal_format(fmt):
    pattern = Abund(pattern="grevesse2007", monh=0)(raw=True)
    encoded = Abund.totype(pattern, fmt, raw=True)
    abundance = Abund(monh=-2, pattern=encoded, type=fmt)
    restored = Abund._load(abundance._save())
    np.testing.assert_allclose(restored(raw=True), abundance(raw=True))
    assert restored.type == fmt
    cloned = copy(abundance)
    np.testing.assert_allclose(cloned(raw=True), abundance(raw=True))
    assert cloned.type == fmt
    cloned.set_A("He", 10.0)
    assert abundance()["He"] == pytest.approx(10.93)
    atmosphere = Atmosphere(monh=-2, abund=encoded, abund_format=fmt)
    extension = atmosphere._save()
    assert extension.header["abund_format"] == "H=12"
    # Also repair previously saved files carrying the misleading input format.
    extension.header["abund_format"] = fmt
    recovered = Atmosphere._load(extension)
    np.testing.assert_allclose(recovered.abund(raw=True), atmosphere.abund(raw=True))


def test_sme_file_roundtrip_preserves_nondefault_abundance_format(tmp_path):
    from pysme.sme import SME_Structure
    sme = SME_Structure()
    grid = make_grid("sme_legacy")
    sme.atmo = grid[0]
    sme.abund = Abund(monh=-2, pattern=legacy_pattern(), type="sme_legacy")
    path = str(tmp_path / "roundtrip.sme")
    sme.save(path)
    restored = SME_Structure.load(path)
    np.testing.assert_allclose(restored.abund(raw=True), sme.abund(raw=True))
    np.testing.assert_allclose(restored.atmo.abund(raw=True), sme.atmo.abund(raw=True))
