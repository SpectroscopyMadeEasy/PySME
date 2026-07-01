# -*- coding: utf-8 -*-
"""Minimal numpy-based kurucz-a1 runtime bundled with PySME."""

from __future__ import annotations

from math import erf, sqrt

import numpy as np

_MODEL_DEPTH_POINTS = 80
_OUTPUT_NAMES = ("RHOX", "T", "P", "XNE", "ABROSS")


def _as_float32_array(value):
    array = np.asarray(value)
    if array.shape == ():
        return np.float32(array.item())
    return array.astype(np.float32, copy=False)


def _unpack_norm_entry(entry):
    if isinstance(entry, np.ndarray) and entry.shape == ():
        entry = entry.item()
    if not isinstance(entry, dict):
        raise TypeError(f"Unexpected normalization entry type {type(entry)!r}")
    return {key: _as_float32_array(value) for key, value in entry.items()}


def load_norm_params(norm_params_path):
    with np.load(norm_params_path, allow_pickle=True) as data:
        return {key: _unpack_norm_entry(value) for key, value in data.items()}


class NormalizationHelper:
    def __init__(self, norm_params):
        self.norm_params = norm_params

    def normalize(self, param_name, data):
        params = self.norm_params[param_name]
        values = np.asarray(data, dtype=np.float32)
        if bool(params["log_scale"]):
            values = np.log10(values + 1e-30)
        return (
            2.0 * (values - params["min"]) / (params["max"] - params["min"]) - 1.0
        ).astype(np.float32, copy=False)

    def denormalize(self, param_name, normalized_data):
        params = self.norm_params[param_name]
        values = np.asarray(normalized_data, dtype=np.float32)
        values = (values + 1.0) / 2.0 * (params["max"] - params["min"]) + params["min"]
        if bool(params["log_scale"]):
            values = np.power(10.0, values) - 1e-30
        return values.astype(np.float32, copy=False)


def _linear(inputs, weight, bias):
    return inputs @ weight.T + bias


def _layer_norm(inputs, weight, bias, eps=1e-5):
    mean = inputs.mean(axis=-1, keepdims=True)
    var = ((inputs - mean) ** 2).mean(axis=-1, keepdims=True)
    return ((inputs - mean) / np.sqrt(var + eps)) * weight + bias


def _gelu(inputs):
    erf_vec = np.vectorize(erf, otypes=[np.float32])
    return 0.5 * inputs * (1.0 + erf_vec(inputs / sqrt(2.0)))


class AtmosphereEmulator:
    def __init__(self, weights, normalizer, default_tau_grid=None):
        self.weights = weights
        self.normalizer = normalizer
        if default_tau_grid is None:
            default_tau_grid = np.logspace(-6, 2, _MODEL_DEPTH_POINTS, dtype=np.float32)
        self.default_tau_grid = np.asarray(default_tau_grid, dtype=np.float32)

    def predict(self, stellar_params, tau_grid=None):
        stellar_params = np.asarray(stellar_params, dtype=np.float32)
        if stellar_params.ndim == 1:
            stellar_params = stellar_params[None, :]

        teff = stellar_params[:, 0:1]
        logg = stellar_params[:, 1:2]
        feh = stellar_params[:, 2:3]
        afe = stellar_params[:, 3:4]

        tau_grid = self._prepare_tau_grid(tau_grid, stellar_params.shape[0])
        current_depth_points = tau_grid.shape[1]
        original_tau_grid = tau_grid.copy()
        if current_depth_points != _MODEL_DEPTH_POINTS:
            tau_grid = self._reshape_tau_grid(tau_grid)

        params_normalized = np.concatenate(
            [
                self.normalizer.normalize("teff", teff),
                self.normalizer.normalize("gravity", logg),
                self.normalizer.normalize("feh", feh),
                self.normalizer.normalize("afe", afe),
                self.normalizer.normalize("TAU", tau_grid),
            ],
            axis=1,
        ).astype(np.float32, copy=False)

        predictions = self._forward(params_normalized)
        if current_depth_points != _MODEL_DEPTH_POINTS:
            predictions = predictions[:, :current_depth_points, :]

        output_features = {
            name: self.normalizer.denormalize(name, predictions[:, :, index])
            for index, name in enumerate(_OUTPUT_NAMES)
        }
        output_features["TAU"] = original_tau_grid

        output_features = {key: value[0] for key, value in output_features.items()}
        output_features["teff"] = float(teff[0, 0])
        output_features["logg"] = float(logg[0, 0])
        output_features["feh"] = float(feh[0, 0])
        output_features["afe"] = float(afe[0, 0])
        output_features["vturb"] = 2.0
        output_features["lonh"] = 1.25
        output_features["geom"] = "PP"
        output_features["citation_info"] = r"""
            @ARTICLE{2025arXiv250706357L,
            author = {{Li}, Jiadong and {Jian}, Mingjie and {Ting}, Yuan-Sen and {Green}, Gregory M.},
            title = "{Differentiable Stellar Atmospheres with Physics-Informed Neural Networks}",
            journal = {arXiv e-prints},
            year = 2025,
            month = jul,
            eid = {arXiv:2507.06357},
            pages = {arXiv:2507.06357},
            doi = {10.48550/arXiv.2507.06357},
            archivePrefix = {arXiv},
            eprint = {2507.06357},
            primaryClass = {astro-ph.SR},
            adsurl = {https://ui.adsabs.harvard.edu/abs/2025arXiv250706357L},
            adsnote = {Provided by the SAO/NASA Astrophysics Data System}
            }
        """
        return output_features

    def _prepare_tau_grid(self, tau_grid, batch_size):
        if tau_grid is None:
            return np.repeat(self.default_tau_grid[None, :], batch_size, axis=0)

        tau_grid = np.asarray(tau_grid, dtype=np.float32)
        if tau_grid.ndim == 1:
            tau_grid = tau_grid[None, :]
        if tau_grid.shape[0] == 1 and batch_size > 1:
            tau_grid = np.repeat(tau_grid, batch_size, axis=0)
        return tau_grid

    @staticmethod
    def _reshape_tau_grid(tau_grid):
        current_depth_points = tau_grid.shape[1]
        if current_depth_points > _MODEL_DEPTH_POINTS:
            return tau_grid[:, :_MODEL_DEPTH_POINTS]

        last_tau = np.repeat(
            tau_grid[:, -1:], _MODEL_DEPTH_POINTS - current_depth_points, axis=1
        )
        return np.concatenate([tau_grid, last_tau], axis=1)

    def _forward(self, inputs):
        batch_size = inputs.shape[0]
        stellar_params = inputs[:, :4]
        tau_values = inputs[:, 4:].reshape(batch_size, _MODEL_DEPTH_POINTS)

        stellar_embedding = self._stellar_encoder(stellar_params)
        tau_embedding = self._tau_encoder(tau_values)
        stellar_embedding = np.repeat(
            stellar_embedding[:, None, :], _MODEL_DEPTH_POINTS, axis=1
        )
        combined = np.concatenate([stellar_embedding, tau_embedding], axis=2)
        return self._predictor(combined)

    def _stellar_encoder(self, values):
        values = _linear(
            values,
            self.weights["stellar_encoder.encoder.0.weight"],
            self.weights["stellar_encoder.encoder.0.bias"],
        )
        values = _layer_norm(
            values,
            self.weights["stellar_encoder.encoder.1.weight"],
            self.weights["stellar_encoder.encoder.1.bias"],
        )
        values = _gelu(values)
        values = _linear(
            values,
            self.weights["stellar_encoder.encoder.3.weight"],
            self.weights["stellar_encoder.encoder.3.bias"],
        )
        values = _layer_norm(
            values,
            self.weights["stellar_encoder.encoder.4.weight"],
            self.weights["stellar_encoder.encoder.4.bias"],
        )
        return _gelu(values)

    def _tau_encoder(self, values):
        values = values.reshape(values.shape[0] * _MODEL_DEPTH_POINTS, 1)
        values = _linear(
            values,
            self.weights["tau_encoder.encoder.0.weight"],
            self.weights["tau_encoder.encoder.0.bias"],
        )
        values = _gelu(values)
        values = _linear(
            values,
            self.weights["tau_encoder.encoder.2.weight"],
            self.weights["tau_encoder.encoder.2.bias"],
        )
        values = _gelu(values)
        return values.reshape(-1, _MODEL_DEPTH_POINTS, values.shape[-1])

    def _predictor(self, values):
        values = _linear(
            values,
            self.weights["predictor.0.weight"],
            self.weights["predictor.0.bias"],
        )
        values = _gelu(values)
        values = _linear(
            values,
            self.weights["predictor.3.weight"],
            self.weights["predictor.3.bias"],
        )
        values = _gelu(values)
        return _linear(
            values,
            self.weights["predictor.6.weight"],
            self.weights["predictor.6.bias"],
        )


def load_from_checkpoint(checkpoint_path, norm_params_path, hidden_size=512, device="cpu"):
    del hidden_size, device
    norm_params = load_norm_params(norm_params_path)
    normalizer = NormalizationHelper(norm_params)
    with np.load(checkpoint_path) as checkpoint:
        weights = {
            key: np.asarray(value, dtype=np.float32) for key, value in checkpoint.items()
        }
    return AtmosphereEmulator(weights, normalizer)
