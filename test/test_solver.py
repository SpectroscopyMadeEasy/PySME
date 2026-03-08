# -*- coding: utf-8 -*-
from os.path import dirname

import numpy as np
import pytest

from pysme.sme import SME_Structure as SME_Struct
from pysme.solve import solve

cwd = dirname(__file__)
filename = "{}/testcase1.inp".format(cwd)


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
