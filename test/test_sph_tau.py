# -*- coding: utf-8 -*-
"""Regression tests for the explicit SPH + TAU model-input path.

The wrapper tests deliberately use a small fake native module.  They check the
ABI selection rules without requiring a locally compiled SMElib, while the
native tests below exercise the actual extension when the new library is
available.
"""

from pathlib import Path

import numpy as np
import pytest

from pysme.sme_synth import SME_DLL


def _atmo(geom="PP", depth_scale="RHOX", ndepth=5):
    """Build a compact, finite atmosphere accepted by InputModelEx."""

    # The native interface expects the input coordinate in increasing order;
    # both RHOX and TAU are represented in that order here.
    depth = np.geomspace(1.0e-6, 1.0e2, ndepth)
    atmosphere = {
        "depth": depth_scale,
        depth_scale: depth,
        "temp": np.linspace(5200.0, 6200.0, ndepth),
        "xne": np.full(ndepth, 1.0e12),
        "xna": np.full(ndepth, 1.0e17),
        "rho": np.full(ndepth, 1.0e-7),
        "vt": np.full(ndepth, 2.0),
        "wlstd": 5000.0,
        "opflag": np.zeros(20, dtype=np.int16),
        "geom": geom,
    }
    if geom == "SPH":
        # Radius is the reference radius and height is the per-layer offset.
        atmosphere["radius"] = 1.0e11
        atmosphere["height"] = np.linspace(5.0e9, 1.0e8, ndepth)
    return atmosphere


class _FakeNative:
    """Capture model calls while exposing the native capability API."""

    def __init__(self, capabilities):
        self.capabilities = capabilities
        self.calls = []

    def GetModelInputCapabilities(self):
        return self.capabilities

    def InputModel(self, *args, **kwargs):
        self.calls.append(("legacy", args, kwargs))

    def InputModelEx(self, *args, **kwargs):
        self.calls.append(("extended", args, kwargs))


def _uninitialised_dll(monkeypatch, native):
    # SME_DLL.__init__ loads data files and the native library.  InputModel is
    # self-contained, so bypass that setup for a focused wrapper test.
    import pysme.sme_synth as sme_synth

    monkeypatch.setattr(sme_synth, "_smelib", native)
    return object.__new__(SME_DLL)


@pytest.mark.parametrize(
    "geom,depth_scale",
    [("PP", "RHOX"), ("PP", "TAU"), ("SPH", "RHOX")],
)
def test_wrapper_keeps_three_legacy_model_inputs(monkeypatch, geom, depth_scale):
    native = _FakeNative("PP:RHOX,PP:TAU,SPH:RHOX")
    dll = _uninitialised_dll(monkeypatch, native)

    dll.InputModel(5750.0, 4.5, 2.0, _atmo(geom, depth_scale))

    assert len(native.calls) == 1
    kind, args, kwargs = native.calls[0]
    assert kind == "legacy"
    assert args[0:4] == (5750.0, 4.5, 5000.0, depth_scale if geom == "PP" else "SPH")
    if geom == "SPH":
        assert "radius" in kwargs
        assert "height" in kwargs
    else:
        assert kwargs == {}


def test_wrapper_dispatches_only_sph_tau_to_extended_abi(monkeypatch):
    native = _FakeNative("PP:RHOX,PP:TAU,SPH:RHOX,SPH:TAU")
    dll = _uninitialised_dll(monkeypatch, native)
    atmo = _atmo("SPH", "TAU")

    dll.InputModel(5750.0, 2.5, 2.0, atmo)

    assert len(native.calls) == 1
    kind, args, kwargs = native.calls[0]
    assert kind == "extended"
    assert args[:5] == (5750.0, 2.5, 5000.0, "SPH", "TAU")
    assert np.array_equal(args[5], atmo["opflag"])
    assert np.array_equal(args[6], atmo["TAU"])
    assert kwargs == {"radius": atmo["radius"], "height": atmo["height"]}


def test_wrapper_rejects_sph_tau_when_old_native_lacks_capability(monkeypatch):
    native = _FakeNative("PP:RHOX,PP:TAU,SPH:RHOX")
    dll = _uninitialised_dll(monkeypatch, native)

    with pytest.raises(RuntimeError, match="SPH:TAU|InputModelEx"):
        dll.InputModel(5750.0, 2.5, 2.0, _atmo("SPH", "TAU"))

    assert native.calls == []


def test_wrapper_rejects_sph_tau_when_extended_symbol_is_missing(monkeypatch):
    class LegacyOnly(_FakeNative):
        InputModelEx = None

    native = LegacyOnly("PP:RHOX,PP:TAU,SPH:RHOX")
    dll = _uninitialised_dll(monkeypatch, native)

    with pytest.raises(RuntimeError, match="SPH:TAU|InputModelEx"):
        dll.InputModel(5750.0, 2.5, 2.0, _atmo("SPH", "TAU"))

    assert native.calls == []


def _native_module_or_skip():
    from pysme.smelib.libtools import get_full_libfile
    from pysme.sme_synth import ensure_smelib_ready

    if not Path(get_full_libfile()).exists():
        pytest.skip("compiled SMElib is not available")
    native = ensure_smelib_ready()
    if not hasattr(native, "InputModelEx"):
        pytest.skip("SMElib extension has no InputModelEx")
    return native


@pytest.mark.parametrize(
    "geom,depth_scale",
    [("PP", "RHOX"), ("PP", "TAU"), ("SPH", "RHOX"), ("SPH", "TAU")],
)
def test_native_accepts_all_explicit_geometry_depth_combinations(geom, depth_scale):
    native = _native_module_or_skip()
    atmo = _atmo(geom, depth_scale)

    native.InputModelEx(
        5750.0,
        2.5 if geom == "SPH" else 4.5,
        atmo["wlstd"],
        geom,
        depth_scale,
        atmo["opflag"],
        atmo[depth_scale],
        atmo["temp"],
        atmo["xne"],
        atmo["xna"],
        atmo["rho"],
        atmo["vt"],
        radius=atmo.get("radius", 0.0),
        height=atmo.get("height"),
    )


def test_native_advertises_the_four_explicit_combinations():
    native = _native_module_or_skip()
    assert native.GetModelInputCapabilities() == (
        "PP:RHOX,PP:TAU,SPH:RHOX,SPH:TAU"
    )


def test_native_can_switch_geometry_and_depth_in_one_process():
    """A mode switch must not retain the previous model's geometry state."""

    native = _native_module_or_skip()
    for geom, depth_scale in (
        ("SPH", "TAU"),
        ("PP", "RHOX"),
        ("SPH", "RHOX"),
        ("PP", "TAU"),
    ):
        atmo = _atmo(geom, depth_scale)
        native.InputModelEx(
            5750.0,
            2.5 if geom == "SPH" else 4.5,
            atmo["wlstd"],
            geom,
            depth_scale,
            atmo["opflag"],
            atmo[depth_scale],
            atmo["temp"],
            atmo["xne"],
            atmo["xna"],
            atmo["rho"],
            atmo["vt"],
            radius=atmo.get("radius", 0.0),
            height=atmo.get("height"),
        )


def test_native_rejects_nonmonotone_explicit_depth():
    native = _native_module_or_skip()
    atmo = _atmo("SPH", "TAU")
    atmo["TAU"] = atmo["TAU"].copy()
    atmo["TAU"][2] = atmo["TAU"][1]

    with pytest.raises((ValueError, RuntimeError), match="depth|monotone|increas"):
        native.InputModelEx(
            5750.0,
            2.5,
            atmo["wlstd"],
            "SPH",
            "TAU",
            atmo["opflag"],
            atmo["TAU"],
            atmo["temp"],
            atmo["xne"],
            atmo["xna"],
            atmo["rho"],
            atmo["vt"],
            radius=atmo["radius"],
            height=atmo["height"],
        )
