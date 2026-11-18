"""Regression tests for the grid/absolute-abundance atmosphere guard."""

import numpy as np
import pytest

from pysme.abund import Abund
from pysme.atmosphere.atmosphere import Atmosphere, AtmosphereError
from pysme.sme import SME_Structure
from pysme.synthesize import Synthesizer


class _RecordingInterpolator:
    def __init__(self):
        self.depth = None
        self.interp = None
        self.geom = None
        self.calls = []

    def interp_atmo_grid(self, source, teff, logg, monh):
        self.calls.append((source, teff, logg, monh))
        atmo = Atmosphere(
            source=source,
            method="grid",
            depth="RHOX",
            interp="TAU",
            geom="PP",
            teff=teff,
            logg=logg,
            monh=monh,
        )
        atmo.interpolation_monh = float(monh)
        atmo.temp = np.array([3200.0, 4800.0])
        return atmo


def _prepared_sme_and_synthesizer():
    sme = SME_Structure(teff=4600, logg=1.1, monh=-3.0)
    sme.atmo.source = "test-grid"
    sme.atmo.method = "grid"
    synthesizer = Synthesizer.__new__(Synthesizer)
    interpolator = _RecordingInterpolator()
    synthesizer.atmosphere_interpolator = interpolator
    synthesizer.get_atmosphere(sme)
    return sme, synthesizer, interpolator


def _atmosphere_structure_bytes(atmo):
    """Stable byte representation of physical fields, excluding method."""
    chunks = []
    for name in atmo.names:
        if name == "method":
            continue
        value = getattr(atmo, name)
        if isinstance(value, np.ndarray):
            chunks.append(value.dtype.str.encode() + value.tobytes())
        else:
            chunks.append(repr(value).encode())
    return b"\0".join(chunks)


def test_absolute_abundance_does_not_trigger_second_grid_interpolation():
    """A frozen metal-poor structure must not be replaced by a solar one."""

    sme, synthesizer, interpolator = _prepared_sme_and_synthesizer()
    first_structure = _atmosphere_structure_bytes(sme.atmo)
    assert interpolator.calls == [("test-grid", 4600, 1.1, -3.0)]
    assert sme.atmo.interpolation_monh == -3.0

    # This is the production failure mode: the absolute pattern's monh is zero
    # even though the atmosphere has already been interpolated at -3.
    sme.abund = Abund(pattern="asplund2009", monh=0.0)
    assert sme.atmo.method == "embedded"
    synthesizer.get_atmosphere(sme)

    assert interpolator.calls == [("test-grid", 4600, 1.1, -3.0)]
    assert _atmosphere_structure_bytes(sme.atmo) == first_structure


def test_in_place_abundance_monh_change_is_rejected():
    """Direct mutation remains guarded because it bypasses ``_abund``."""

    sme, synthesizer, interpolator = _prepared_sme_and_synthesizer()
    sme.abund.monh = 0.0
    with pytest.raises(AtmosphereError, match="already interpolated"):
        synthesizer.get_atmosphere(sme)
    assert sme.atmo.method == "grid"
    assert len(interpolator.calls) == 1


def test_public_monh_change_requests_a_new_grid_interpolation():
    """The ordinary ``sme.monh = value`` workflow keeps its old behavior."""

    sme, synthesizer, interpolator = _prepared_sme_and_synthesizer()
    sme.monh = -2.0
    assert sme.atmo.method == "grid"
    assert sme.atmo.interpolation_monh is None

    synthesizer.get_atmosphere(sme)
    assert interpolator.calls == [
        ("test-grid", 4600, 1.1, -3.0),
        ("test-grid", 4600, 1.1, -2.0),
    ]
    assert sme.atmo.interpolation_monh == -2.0


def test_consistent_abundance_replacement_keeps_grid_behavior():
    """A replacement at the recorded coordinate follows the existing path."""

    sme, synthesizer, interpolator = _prepared_sme_and_synthesizer()
    first_structure = _atmosphere_structure_bytes(sme.atmo)
    sme.abund = Abund(pattern="asplund2009", monh=-3.0)
    assert sme.atmo.method == "grid"

    synthesizer.get_atmosphere(sme)
    assert interpolator.calls == [
        ("test-grid", 4600, 1.1, -3.0),
        ("test-grid", 4600, 1.1, -3.0),
    ]
    assert _atmosphere_structure_bytes(sme.atmo) == first_structure
