# -*- coding: utf-8 -*-
import numpy as np

from pysme.abund import Abund
from pysme.atmosphere.atmosphere import Atmosphere
from pysme.atmosphere.kurucz_a1 import load_from_checkpoint, load_norm_params
from pysme.atmosphere.providers import KuruczA1Provider, resolve_routine_atmosphere_provider
from pysme.sme import SME_Structure
from pysme.synthesize import Synthesizer


class _FakeEmulator:
    def __init__(self):
        self.calls = []

    def predict(self, stellar_params, tau_grid=None):
        self.calls.append((stellar_params, tau_grid))
        return {
            "RHOX": np.array([1e-5, 3e-5, 1e-4], dtype=float),
            "TAU": np.array([1e-6, 3e-6, 1e-5], dtype=float),
            "T": np.array([4500.0, 4700.0, 5000.0], dtype=float),
            "P": np.array([1e2, 2e2, 4e2], dtype=float),
            "XNE": np.array([1e10, 2e10, 4e10], dtype=float),
            "teff": 5000.0,
            "logg": 4.4,
            "vturb": 2.0,
            "lonh": 1.25,
            "geom": "PP",
            "citation_info": "kurucz-a1-citation",
        }


def test_kurucz_a1_provider_converts_prediction():
    provider = KuruczA1Provider()
    provider._emulator = _FakeEmulator()

    sme = SME_Structure()
    sme.teff = 5000.0
    sme.logg = 4.4
    sme.monh = -0.2
    sme.abund = Abund(monh=-0.2, pattern="asplund2009")
    sme.abund.set_xm("Mg", 0.3)
    sme.abund.set_xm("Si", 0.3)
    atmo = Atmosphere(method="routine", source="kurucz-a1", depth="RHOX", interp="TAU")

    out = provider(sme, atmo)

    assert np.allclose(out.rhox, [1e-5, 3e-5, 1e-4])
    assert np.allclose(out.tau, [1e-6, 3e-6, 1e-5])
    assert np.allclose(out.temp, [4500.0, 4700.0, 5000.0])
    assert np.allclose(out.xne, [1e10, 2e10, 4e10])
    assert np.all(np.isfinite(out.xna))
    assert np.all(np.isfinite(out.rho))
    assert out.method == "embedded"
    assert out.source == "kurucz-a1"
    assert out.citation_info == "kurucz-a1-citation"


def test_kurucz_a1_provider_bounds_from_norm_params():
    provider = KuruczA1Provider(lfs_atmo=object())
    provider._norm_params = {
        "teff": {"min": np.log10(3500.0), "max": np.log10(8000.0), "log_scale": True},
        "gravity": {"min": -0.5, "max": 5.0, "log_scale": False},
        "feh": {"min": -2.5, "max": 0.5, "log_scale": False},
    }

    bounds = provider.get_bounds()

    assert np.isclose(bounds["teff"][0], 3500.0)
    assert np.isclose(bounds["teff"][1], 8000.0)
    assert bounds["logg"] == (-0.5, 5.0)
    assert bounds["monh"] == (-2.5, 0.5)


def test_synthesizer_resolves_named_routine_provider(monkeypatch):
    fake_provider = KuruczA1Provider(lfs_atmo=object())
    fake_provider._emulator = _FakeEmulator()

    monkeypatch.setattr(
        "pysme.synthesize.resolve_routine_atmosphere_provider",
        lambda source: fake_provider,
    )

    sme = SME_Structure()
    sme.teff = 5000.0
    sme.logg = 4.4
    sme.abund = Abund(monh=0.0, pattern="asplund2009")
    sme.atmo.method = "routine"
    sme.atmo.source = "kurucz-a1"

    synth = Synthesizer()
    out = synth.get_atmosphere(sme)

    assert out.atmo.method == "embedded"
    assert np.allclose(out.atmo.temp, [4500.0, 4700.0, 5000.0])


def test_resolve_named_provider_returns_kurucz_a1():
    provider = resolve_routine_atmosphere_provider("kurucz-a1")
    assert isinstance(provider, KuruczA1Provider)


def test_kurucz_a1_provider_uses_lfs_paths():
    class _FakeLFS:
        def __init__(self):
            self.requests = []

        def get(self, key):
            self.requests.append(key)
            return f"/tmp/{key}"

    provider = KuruczA1Provider(lfs_atmo=_FakeLFS())

    assert provider._resolve_norm_params_path() == "/tmp/kurucz_a1_norm_params.npz"
    assert provider._resolve_checkpoint_path() == "/tmp/kurucz_a1_weights.npz"


def test_kurucz_a1_numpy_loader_round_trip(tmp_path):
    norm_path = tmp_path / "norm.npz"
    weights_path = tmp_path / "weights.npz"
    embed = 2
    combined = embed * 2

    np.savez(
        norm_path,
        teff=np.array(
            {"min": np.float32(3.0), "max": np.float32(4.0), "log_scale": True},
            dtype=object,
        ),
        gravity=np.array(
            {"min": np.float32(0.0), "max": np.float32(5.0), "log_scale": False},
            dtype=object,
        ),
        feh=np.array(
            {"min": np.float32(-2.0), "max": np.float32(0.5), "log_scale": False},
            dtype=object,
        ),
        afe=np.array(
            {"min": np.float32(-0.5), "max": np.float32(0.5), "log_scale": False},
            dtype=object,
        ),
        TAU=np.array(
            {"min": np.float32(-6.0), "max": np.float32(2.0), "log_scale": True},
            dtype=object,
        ),
        RHOX=np.array(
            {"min": np.float32(-10.0), "max": np.float32(-2.0), "log_scale": True},
            dtype=object,
        ),
        T=np.array(
            {"min": np.float32(3.0), "max": np.float32(4.5), "log_scale": True},
            dtype=object,
        ),
        P=np.array(
            {"min": np.float32(0.0), "max": np.float32(8.0), "log_scale": True},
            dtype=object,
        ),
        XNE=np.array(
            {"min": np.float32(8.0), "max": np.float32(16.0), "log_scale": True},
            dtype=object,
        ),
        ABROSS=np.array(
            {"min": np.float32(-6.0), "max": np.float32(2.0), "log_scale": True},
            dtype=object,
        ),
    )
    np.savez(
        weights_path,
        **{
            "stellar_encoder.encoder.0.weight": np.zeros((embed, 4), dtype=np.float32),
            "stellar_encoder.encoder.0.bias": np.zeros(embed, dtype=np.float32),
            "stellar_encoder.encoder.1.weight": np.ones(embed, dtype=np.float32),
            "stellar_encoder.encoder.1.bias": np.zeros(embed, dtype=np.float32),
            "stellar_encoder.encoder.3.weight": np.zeros((embed, embed), dtype=np.float32),
            "stellar_encoder.encoder.3.bias": np.zeros(embed, dtype=np.float32),
            "stellar_encoder.encoder.4.weight": np.ones(embed, dtype=np.float32),
            "stellar_encoder.encoder.4.bias": np.zeros(embed, dtype=np.float32),
            "tau_encoder.encoder.0.weight": np.zeros((embed // 2, 1), dtype=np.float32),
            "tau_encoder.encoder.0.bias": np.zeros(embed // 2, dtype=np.float32),
            "tau_encoder.encoder.2.weight": np.zeros((embed, embed // 2), dtype=np.float32),
            "tau_encoder.encoder.2.bias": np.zeros(embed, dtype=np.float32),
            "predictor.0.weight": np.zeros((combined * 2, combined), dtype=np.float32),
            "predictor.0.bias": np.zeros(combined * 2, dtype=np.float32),
            "predictor.3.weight": np.zeros((combined * 2, combined * 2), dtype=np.float32),
            "predictor.3.bias": np.zeros(combined * 2, dtype=np.float32),
            "predictor.6.weight": np.zeros((5, combined * 2), dtype=np.float32),
            "predictor.6.bias": np.zeros(5, dtype=np.float32),
        },
    )

    norm = load_norm_params(norm_path)
    emulator = load_from_checkpoint(weights_path, norm_path)
    out = emulator.predict([5000.0, 4.0, 0.0, 0.0], tau_grid=np.logspace(-5, 1, 7))

    assert bool(norm["teff"]["log_scale"]) is True
    assert out["TAU"].shape == (7,)
    assert out["RHOX"].shape == (7,)
    assert np.all(np.isfinite(out["T"]))
