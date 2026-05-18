# -*- coding: utf-8 -*-
from os.path import dirname

import numpy as np
import pytest

from pysme.abund import Abund
from pysme.sme import SME_Structure as SME_Struct
from pysme.solve import SME_Solver, solve

cwd = dirname(__file__)
filename = "{}/testcase1.inp".format(cwd)


class _DummyProgressBar:
    def __init__(self):
        self.total = 0

    def update(self, *_args, **_kwargs):
        return None


def _make_free_abundance_sme():
    sme = SME_Struct()
    sme.monh = -1.196
    sme.abund = Abund(monh=sme.monh, pattern="asplund2021")

    effective_ti = 3.909
    pattern_ti = effective_ti - sme.monh
    sme["Abund Ti"] = pattern_ti

    assert sme.abund["Ti"] == pytest.approx(effective_ti)
    assert sme.abund.get_pattern(type="H=12")["Ti"] == pytest.approx(pattern_ti)

    return sme, pattern_ti, effective_ti


def test_simple():
    sme = SME_Struct.load(filename)
    sme2 = solve(sme, ["teff"])

    assert sme2.synth is not None
    assert sme2.fitresults is not None
    assert sme2.fitresults.covariance is not None
    assert isinstance(sme2.fitresults.covariance, np.ndarray)
    assert np.all(sme2.fitresults.covariance != 0)

    assert isinstance(sme2.fitresults.uncertainties, np.ndarray)
    assert len(sme2.fitresults.uncertainties) == 1
    assert sme2.fitresults.parameters[0] == "teff"
    assert sme2.fitresults.uncertainties[0] != 0

    assert np.array_equal(sme2.fitresults.covariance.shape, [1, 1])
    assert sme2.fitresults.covariance.ndim == 2

    assert sme2.fitresults.chisq is not None
    assert sme2.fitresults.chisq != 0


def test_solve_requires_wave():
    sme = SME_Struct()
    sme.spec = [np.ones(10)]
    with pytest.raises(ValueError, match="wavelength grid"):
        solve(sme, ["teff"])


def test_solve_requires_spec():
    sme = SME_Struct()
    sme.wave = [np.linspace(5000.0, 5001.0, 10)]
    with pytest.raises(ValueError, match="observed spectrum"):
        solve(sme, ["teff"])


@pytest.mark.parametrize("field_name", ["uncs", "mask"])
def test_solve_rejects_mismatched_segment_length(field_name):
    sme = SME_Struct.load(filename)
    bad = np.ones(max(1, len(sme.wave[0]) - 1))
    if field_name == "mask":
        bad = bad.astype(int)
    setattr(sme, field_name, [bad])

    with pytest.raises(ValueError, match=field_name):
        solve(sme, ["teff"])


def test_solver_default_free_abundance_uses_pattern_scale():
    sme, pattern_ti, effective_ti = _make_free_abundance_sme()
    solver = SME_Solver()
    solver.parameter_names = ["Abund Ti"]

    p0 = solver.get_default_values(sme)

    assert p0[0] == pytest.approx(pattern_ti)
    assert p0[0] != pytest.approx(effective_ti)


def test_solver_residuals_keep_effective_abundance_at_initial_guess(monkeypatch):
    sme, _pattern_ti, effective_ti = _make_free_abundance_sme()
    solver = SME_Solver()
    solver.parameter_names = ["Abund Ti"]
    solver.progressbar = _DummyProgressBar()
    solver.progressbar_jacobian = _DummyProgressBar()

    captured = {}

    def fake_synthesize_spectrum(local_sme, **_kwargs):
        captured["effective_ti"] = local_sme.abund["Ti"]
        raise RuntimeError("stop after assignment")

    monkeypatch.setattr(solver.synthesizer, "synthesize_spectrum", fake_synthesize_spectrum)

    p0 = solver.get_default_values(sme)
    with pytest.raises(RuntimeError, match="stop after assignment"):
        solver._residuals(
            p0,
            sme,
            spec=np.array([[0.0]]),
            uncs=np.array([[1.0]]),
            mask=np.array([True]),
            segments=[0],
        )

    assert captured["effective_ti"] == pytest.approx(effective_ti)


def test_solver_backup_restore_preserves_effective_abundance(tmp_path):
    sme, _pattern_ti, effective_ti = _make_free_abundance_sme()
    solver = SME_Solver(filename=str(tmp_path / "solver_state.sme"), restore=True)
    solver.parameter_names = ["Abund Ti"]

    solver.backup(sme)

    sme["Abund Ti"] = 0.0
    assert sme.abund["Ti"] != pytest.approx(effective_ti)

    restored = solver.restore_func(sme)

    assert restored.abund["Ti"] == pytest.approx(effective_ti)
