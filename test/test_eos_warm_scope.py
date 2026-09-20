import pytest
import numpy as np

import pysme.solve as solve_module


class _WarmControl:
    def __init__(self):
        self.calls = []

    def SetEosWarmStartMode(self, mode):
        self.calls.append(bool(mode))


class _Synthesizer:
    def __init__(self, dll):
        self._dll = dll

    def get_dll(self):
        return self._dll


def _install_fake_solver(monkeypatch, *, fail=False):
    control = _WarmControl()

    class _Solver:
        def __init__(self, filename=None, restore=False):
            self.synthesizer = _Synthesizer(control)

        def solve(self, sme, param_names, segments, **kwargs):
            if fail:
                raise RuntimeError("fit failed")
            return sme

    monkeypatch.setattr(solve_module, "SME_Solver", _Solver)
    return control


def test_abundance_fit_scopes_exact_eos_history(monkeypatch):
    control = _install_fake_solver(monkeypatch)
    sme = object()

    assert solve_module.solve(sme, ["Abund Fe"]) is sme
    assert control.calls == [True, False]


def test_non_eos_fit_does_not_enable_eos_history(monkeypatch):
    control = _install_fake_solver(monkeypatch)

    solve_module.solve(object(), ["vrad"])
    assert control.calls == []


def test_fitparameters_array_scopes_exact_eos_history(monkeypatch):
    control = _install_fake_solver(monkeypatch)

    class _Sme:
        fitparameters = np.array(["Abund Fe"])

    sme = _Sme()
    assert solve_module.solve(sme) is sme
    assert control.calls == [True, False]


def test_abundance_fit_clears_eos_history_after_failure(monkeypatch):
    control = _install_fake_solver(monkeypatch, fail=True)

    with pytest.raises(RuntimeError, match="fit failed"):
        solve_module.solve(object(), ["Abund Fe"])
    assert control.calls == [True, False]
