# -*- coding: utf-8 -*-
"""Named dynamic atmosphere providers."""

from __future__ import annotations

import numpy as np

from ..abund import _atom_weight
from ..large_file_storage import setup_atmo
from .atmosphere import Atmosphere, AtmosphereError

_KB = 1.380649e-16
_MH = 1.6735575e-24
_ALPHA_ELEMENTS = ("O", "Mg", "Si", "Ca", "Ti")
_PROVIDER_CACHE = {}


def resolve_routine_atmosphere_provider(source):
    """Resolve a routine-atmosphere provider from a string or callable."""
    if callable(source):
        return source
    if not isinstance(source, str):
        raise AtmosphereError(
            "Routine atmosphere source must be a callable or a known provider name."
        )

    key = source.strip().casefold()
    if key in {"kurucz-a1", "kuruczone"}:
        provider = _PROVIDER_CACHE.get("kurucz-a1")
        if provider is None:
            provider = KuruczA1Provider()
            _PROVIDER_CACHE["kurucz-a1"] = provider
        return provider

    raise AtmosphereError(f"Unknown routine atmosphere provider {source!r}")


class KuruczA1Provider:
    """Adapter that turns kurucz-a1 emulator output into a PySME atmosphere."""

    def __init__(
        self,
        lfs_atmo=None,
        checkpoint_key="kurucz_a1_weights.npz",
        norm_params_key="kurucz_a1_norm_params.npz",
        checkpoint_path=None,
        norm_params_path=None,
    ):
        self.lfs_atmo = setup_atmo() if lfs_atmo is None else lfs_atmo
        self.checkpoint_key = checkpoint_key
        self.norm_params_key = norm_params_key
        self.checkpoint_path = checkpoint_path
        self.norm_params_path = norm_params_path
        self._emulator = None
        self._norm_params = None

    def __call__(self, sme, atmo):
        emulator = self._get_emulator()
        alpha_enhancement = self._get_alpha_enhancement(sme)
        tau_grid = None
        if atmo.tau is not None:
            tau_grid = np.asarray(atmo.tau, dtype=float)

        prediction = emulator.predict(
            [sme.teff, sme.logg, sme.monh, alpha_enhancement], tau_grid=tau_grid
        )
        return self._prediction_to_atmosphere(sme, atmo, prediction)

    def get_bounds(self):
        norm_params = self._get_norm_params()
        return {
            "teff": self._param_bounds(norm_params["teff"]),
            "logg": self._param_bounds(norm_params["gravity"]),
            "monh": self._param_bounds(norm_params["feh"]),
        }

    def _prediction_to_atmosphere(self, sme, atmo, prediction):
        atmosphere = Atmosphere()
        atmosphere.teff = float(prediction.get("teff", sme.teff))
        atmosphere.logg = float(prediction.get("logg", sme.logg))
        atmosphere.abund = sme.abund.__copy__()
        atmosphere.vturb = float(prediction.get("vturb", atmo.vturb))
        atmosphere.lonh = float(prediction.get("lonh", atmo.lonh))
        atmosphere.source = "kurucz-a1"
        atmosphere.method = "embedded"
        atmosphere.geom = prediction.get("geom", atmo.geom if atmo.geom else "PP")
        atmosphere.depth = atmo.depth if atmo.depth is not None else "RHOX"
        atmosphere.interp = atmo.interp if atmo.interp is not None else "TAU"
        atmosphere.wlstd = atmo.wlstd
        atmosphere.opflag = np.array(atmo.opflag, copy=True)
        atmosphere.rhox = np.asarray(prediction["RHOX"], dtype=float)
        atmosphere.tau = np.asarray(prediction["TAU"], dtype=float)
        atmosphere.temp = np.asarray(prediction["T"], dtype=float)
        atmosphere.xne = np.asarray(prediction["XNE"], dtype=float)

        pressure = np.asarray(prediction["P"], dtype=float)
        atmosphere.xna = pressure / (_KB * atmosphere.temp)
        mu = self._mean_particle_weight(sme.abund)
        atmosphere.rho = pressure * mu * _MH / (_KB * atmosphere.temp)

        citation = prediction.get("citation_info", "")
        if citation:
            atmosphere.citation_info = citation
        return atmosphere

    def _get_alpha_enhancement(self, sme):
        values = []
        for element in _ALPHA_ELEMENTS:
            try:
                values.append(float(sme.abund.get_xm(element)))
            except (KeyError, ValueError):
                continue
        if not values:
            return 0.0
        finite = np.asarray(values, dtype=float)
        finite = finite[np.isfinite(finite)]
        if finite.size == 0:
            return 0.0
        return float(np.mean(finite))

    def _mean_particle_weight(self, abund):
        ratios = abund.get_pattern(type="n/nH", raw=True)
        valid = np.isfinite(ratios)
        total_number = float(np.nansum(ratios[valid]))
        if total_number <= 0:
            raise AtmosphereError(
                "Abundance pattern does not define a positive total particle number."
            )
        total_mass = float(np.nansum(ratios[valid] * np.asarray(_atom_weight)[valid]))
        return total_mass / total_number

    def _param_bounds(self, params):
        minimum = self._undo_transform(params["min"], params["log_scale"])
        maximum = self._undo_transform(params["max"], params["log_scale"])
        return float(minimum), float(maximum)

    @staticmethod
    def _undo_transform(value, log_scale):
        if hasattr(value, "item"):
            value = value.item()
        value = float(value)
        return 10**value if log_scale else value

    def _get_norm_params(self):
        if self._norm_params is not None:
            return self._norm_params

        from .kurucz_a1 import load_norm_params

        norm_params_path = self._resolve_norm_params_path()
        self._norm_params = load_norm_params(norm_params_path)
        return self._norm_params

    def _get_emulator(self):
        if self._emulator is not None:
            return self._emulator

        from .kurucz_a1 import load_from_checkpoint

        checkpoint_path = self._resolve_checkpoint_path()
        norm_params_path = self._resolve_norm_params_path()
        self._emulator = load_from_checkpoint(
            checkpoint_path, norm_params_path, device="cpu"
        )
        return self._emulator

    def _resolve_checkpoint_path(self):
        if self.checkpoint_path is not None:
            return self.checkpoint_path
        return self.lfs_atmo.get(self.checkpoint_key)

    def _resolve_norm_params_path(self):
        if self.norm_params_path is not None:
            return self.norm_params_path
        return self.lfs_atmo.get(self.norm_params_key)
