# -*- coding: utf-8 -*-
from collections import OrderedDict
import warnings

import numpy as np
import pytest

from pysme.abund import Abund
from pysme.sme import SME_Structure

pattern_names = ["Asplund2009", "Grevesse2007", "Empty"]
types = ["H=12", "n/nH", "n/nTot", "SME"]


def test_init_with_too_few_args():
    """Test that __init__ raise an error if too few arguments are passed."""
    Abund()
    Abund(monh=0)
    Abund(pattern="asplund2009")
    Abund(type="H=12")


def test_init_using_pattern_names():
    """Test handling of abundance pattern name passed to __init__()."""
    # Each abundance pattern name yields an Abund object.
    for pattern_name in pattern_names:
        assert isinstance(Abund(pattern=pattern_name, monh=0), Abund)

    # The 'Empty' abundance pattern has a value of None for all elements.
    abund = Abund(monh=0, pattern="Empty")
    pattern = abund.get_pattern(raw=True)
    assert pattern[0] == 0
    assert np.all(np.isnan(pattern[1:]))

    # An invalid abundance pattern name raises an error.
    with pytest.raises(ValueError):
        Abund(monh=0, pattern="INVALID")


def test_call_returns_abund_in_odict():
    """Test return value, which is an ordered dictionary with element
    abbreviations as the keys and abundances as the values.
    """
    abund = Abund(pattern=pattern_names[0], monh=0)
    assert isinstance(abund(), dict)
    assert tuple(abund().keys()) == abund.elem


def test_getitem_returns_abund_values():
    """Test getitem method, which return computed abundance values for
    the specified element or list of elements.
    """
    abund = Abund(pattern=pattern_names[0], monh=0)
    assert abund["H"] == 12


def test_monh_property_set_and_get():
    """Test setting and getting monh property. Set converts input to float."""
    # Input str convertable to float yields a float with the specified value.
    abund = Abund(pattern=pattern_names[0], monh="-6e-1")
    assert isinstance(abund.monh, float)
    assert abund.monh == -0.6

    # Input int yields a float with the specified value.
    abund.monh = -2
    assert isinstance(abund.monh, float)
    assert abund.monh == -2.0

    # Input float yields a float with the specified value.
    abund.monh = 0.3
    assert isinstance(abund.monh, float)
    assert abund.monh == 0.3

    # Input str that cannot be converted to float raises an error.
    with pytest.raises(ValueError):
        abund = Abund(pattern=pattern_names[0], monh="ABC")

    # Input that is not a string or a number raises an error.
    with pytest.raises(TypeError):
        abund.monh = []


def test_pattern_property_set_and_get():
    """Test setting and getting pattern property. Set is not allowed."""
    # Raise error is user tries to set pattern
    abund = Abund(pattern="Empty", monh=0)
    with pytest.raises(AttributeError):
        abund.pattern = 0.0
    assert isinstance(dict(abund.pattern), dict)


def test_update_pattern():
    """Test behavior of update_pattern(), which modifies values in _pattern
    for the specified element(s).
    """
    # Update for one element yields float with the specified value.
    abund = Abund(pattern="Empty", monh=0)
    assert np.isnan(abund["Fe"])
    abund.update_pattern({"Fe": "3.14"})
    assert isinstance(abund["Fe"], float)
    assert abund["Fe"] == 3.14

    # Update for two elements yields floats with the specified values.
    abund.update_pattern({"C": 8.4, "F": 5})
    assert isinstance(abund["C"], float)
    assert isinstance(abund["F"], float)
    assert abund["C"] == 8.4
    assert abund["F"] == 5.0


def test_totype_fromtype():
    """Test behavior of totype() and fromtype(), which are static methods
    that return a copy of the input abundance pattern transformed to or
    from the specified abudnance pattern type.
    """
    # Round trip tests that compare copy=fromtype(totype()) with original.
    orig = Abund(pattern=pattern_names[0], monh=0)()
    for type in types:
        pattern = Abund.totype(orig, type)
        copy = Abund.fromtype(pattern, type)
        # Same elements in the same order for full dictionary.
        assert copy.keys() == orig.keys()
        # Same elements have abundance defined (!= None).
        o = OrderedDict((k, v) for k, v in orig.items() if not np.isnan(v))
        c = OrderedDict((k, v) for k, v in copy.items() if not np.isnan(v))
        assert c.keys() == o.keys()
        # Logarithmic abundances differ by less than 1e-10.
        assert all([abs(o[k] - c[k]) < 1e-10 for k in o.keys()])
        # Lowercase type yields same result as mixed case type.
        # type_lc = type.lower()
        # assert copy == Abund.fromtype(Abund.totype(orig, type_lc), type_lc)

    # Invalid abundance pattern type raises error.
    with pytest.raises(ValueError):
        copy = Abund.totype(orig, "INVALID")
    with pytest.raises(ValueError):
        copy = Abund.fromtype(orig, "INVALID")


def _make_asplund2021_abund():
    abund = Abund(pattern="asplund2021", monh=-1.196)
    abund.xm["Ti"] = 0.135
    return abund


def test_abundance_views_basic_consistency():
    abund = _make_asplund2021_abund()

    assert abund.reference["Ti"] == pytest.approx(4.97)
    assert abund.pattern["Ti"] == pytest.approx(5.105)
    assert abund.A["Ti"] == pytest.approx(3.909)
    assert abund.xh["Ti"] == pytest.approx(-1.061)
    assert abund.xm["Ti"] == pytest.approx(0.135)


@pytest.mark.parametrize(
    ("setter", "value"),
    [
        ("A", 3.909),
        ("xh", -1.061),
        ("xm", 0.135),
        ("pattern", 5.105),
    ],
)
def test_abundance_view_setters_are_equivalent(setter, value):
    abund = Abund(pattern="asplund2021", monh=-1.196)
    getattr(abund, setter)["Ti"] = value

    assert abund.A["Ti"] == pytest.approx(3.909)
    assert abund.xm["Ti"] == pytest.approx(0.135)
    assert abund.pattern["Ti"] == pytest.approx(5.105)


def test_changing_monh_preserves_xm_and_pattern():
    abund = _make_asplund2021_abund()
    old_pattern = abund.pattern["Ti"]

    abund.monh = -1.3

    assert abund.pattern["Ti"] == pytest.approx(old_pattern)
    assert abund.xm["Ti"] == pytest.approx(0.135)
    assert abund.A["Ti"] == pytest.approx(old_pattern - 1.3)


def test_reference_remains_fixed_after_updates():
    abund = Abund(pattern="asplund2021", monh=-1.196)
    ref = abund.reference["Ti"]

    abund.A["Ti"] = 3.909
    abund.xh["Ti"] = -1.061
    abund.xm["Ti"] = 0.135
    abund.pattern["Ti"] = 5.105

    assert abund.reference["Ti"] == pytest.approx(ref)


def test_solar_alias_for_builtin_and_error_for_custom_pattern():
    solar = Abund(pattern="asplund2021", monh=0)
    assert solar.solar["Ti"] == pytest.approx(solar.reference["Ti"])

    custom = Abund(pattern=solar.get_pattern("H=12", raw=True), monh=0, type="H=12")
    assert custom.reference["Ti"] == pytest.approx(solar.reference["Ti"])
    with pytest.raises(ValueError, match="not marked as a solar pattern"):
        _ = custom.solar["Ti"]


@pytest.mark.parametrize("elem", ["H", "He"])
def test_xm_invalid_for_h_and_he(elem):
    abund = Abund(pattern="asplund2021", monh=-1.196)
    with pytest.raises(ValueError, match="only defined for elements heavier than He"):
        _ = abund.xm[elem]
    with pytest.raises(ValueError, match="only defined for elements heavier than He"):
        abund.xm[elem] = 0.0


def test_legacy_setter_warns_and_pattern_view_does_not():
    abund = Abund(pattern="asplund2021", monh=-1.196)

    with pytest.warns(FutureWarning, match="Direct abundance assignment"):
        abund["Ti"] = 5.105
    assert abund.pattern["Ti"] == pytest.approx(5.105)

    with warnings.catch_warnings(record=True) as record:
        warnings.simplefilter("always")
        abund.pattern["Ti"] = 5.000
    assert not record
    assert abund.pattern["Ti"] == pytest.approx(5.000)


def test_sme_structure_abundance_assignment_avoids_warning():
    sme = SME_Structure()
    sme.abund = Abund(pattern="asplund2021", monh=-1.196)

    with warnings.catch_warnings(record=True) as record:
        warnings.simplefilter("always")
        sme["Abund Ti"] = 5.105
    assert not record
    assert sme.abund.pattern["Ti"] == pytest.approx(5.105)


def test_no_xfe_view_added():
    abund = Abund(pattern="asplund2021", monh=0)
    assert not hasattr(abund, "xfe")


def test_abundance_reference_metadata_roundtrip():
    abund = _make_asplund2021_abund()

    restored = Abund.from_dict(abund.to_dict())

    assert restored.reference["Ti"] == pytest.approx(abund.reference["Ti"])
    assert restored.pattern["Ti"] == pytest.approx(abund.pattern["Ti"])
    assert restored.A["Ti"] == pytest.approx(abund.A["Ti"])
