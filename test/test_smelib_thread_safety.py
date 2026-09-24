from concurrent.futures import ThreadPoolExecutor
from copy import deepcopy
from pathlib import Path
from threading import Barrier, Event

import numpy as np
import pytest

from pysme.abund import Abund
from pysme.atmosphere.krzfile import KrzFile
from pysme.linelist.linelist import LineList
from pysme.sme_synth import SME_DLL

from .conftest import skipif_smelib

pytestmark = pytest.mark.filterwarnings(
    "ignore:.*EOS-computed electron density differs from the model"
)


def _make_linelist():
    linelist = LineList()
    linelist.add("Fe 1", 5502.9931, 0.9582, -3.047, 7.19, -6.22, 239.249)
    linelist.add("Cr 2", 5503.5955, 4.1682, -2.117, 8.37, -6.49, 195.248)
    return linelist


def _make_cases():
    atmosphere_a = KrzFile(str(Path(__file__).with_name("testatmo1.krz")))
    atmosphere_b = deepcopy(atmosphere_a)
    atmosphere_b.temp = np.asarray(atmosphere_b.temp) * 0.82

    abundance_a = Abund(monh=-0.5, pattern="Asplund2009")
    abundance_b = Abund(monh=0.5, pattern="Asplund2009")
    abundance_b.A["Fe"] += 0.8

    return (
        (5750.0, atmosphere_a, abundance_a),
        (4700.0, atmosphere_b, abundance_b),
    )


def _configure(dll, linelist, case):
    teff, atmosphere, abundance = case
    dll.SetLibraryPath()
    dll.InputLineList(linelist)
    dll.InputModel(teff, 4.5, 2.0, atmosphere)
    dll.InputAbund(abundance)
    dll.Ionization(7)
    dll.SetVWscale(2.5)
    dll.SetH2broad()
    dll.InputWaveRange(5502.5, 5504.2)
    dll.Opacity()


def _observe(dll):
    _, wave, flux, continuum = dll.Transf([1.0], accrt=0.1, accwi=0.1)
    return {
        "wave": np.asarray(wave).copy(),
        "flux": np.asarray(flux[0]).copy(),
        "continuum": np.asarray(continuum[0]).copy(),
        "density": np.asarray(dll.GetDensity()).copy(),
        "natom": np.asarray(dll.GetNatom()).copy(),
        "nelec": np.asarray(dll.GetNelec()).copy(),
    }


def _run_case(dll, linelist, case):
    _configure(dll, linelist, case)
    return _observe(dll)


def _assert_observation_matches(actual, expected):
    assert actual.keys() == expected.keys()
    for name in actual:
        assert actual[name].shape == expected[name].shape, name
        assert np.allclose(actual[name], expected[name], rtol=1e-10, atol=1e-12), name


def _observation_matches(actual, expected, name):
    return actual[name].shape == expected[name].shape and np.allclose(
        actual[name], expected[name], rtol=1e-10, atol=1e-12
    )


def _mismatch_counts(observations, expected):
    return {
        name: sum(
            not _observation_matches(actual, expected, name)
            for actual in observations
        )
        for name in expected
    }


@skipif_smelib
def test_smelib_serial_abab_is_stable():
    linelist = _make_linelist()
    case_a, case_b = _make_cases()
    dll_a = SME_DLL()
    dll_b = SME_DLL()
    dll_a.SetEosWarmStartMode(True)
    try:
        reference_a = _run_case(dll_a, linelist, case_a)
        reference_b = _run_case(dll_b, linelist, case_b)
        repeated_a = _run_case(dll_a, linelist, case_a)
        repeated_b = _run_case(dll_b, linelist, case_b)
    finally:
        dll_a.SetEosWarmStartMode(False)

    _assert_observation_matches(repeated_a, reference_a)
    _assert_observation_matches(repeated_b, reference_b)


@skipif_smelib
@pytest.mark.xfail(
    strict=True,
    reason="SMElib synthesis state is process-global and not isolated by SME_DLL instance",
)
def test_interleaved_smelib_instances_preserve_flux_and_eos_outputs():
    linelist = _make_linelist()
    case_a, case_b = _make_cases()
    dll_a = SME_DLL()
    dll_b = SME_DLL()
    dll_a.SetEosWarmStartMode(True)

    try:
        reference_a = _run_case(dll_a, linelist, case_a)
        reference_b = _run_case(dll_b, linelist, case_b)

        repetitions = 20
        a_configured = [Event() for _ in range(repetitions)]
        b_configured = [Event() for _ in range(repetitions)]
        a_observed = [Event() for _ in range(repetitions)]
        b_observed = [Event() for _ in range(repetitions)]

        def run_a():
            observations = []
            for index in range(repetitions):
                _configure(dll_a, linelist, case_a)
                a_configured[index].set()
                assert b_configured[index].wait(timeout=10)
                observations.append(_observe(dll_a))
                a_observed[index].set()
                assert b_observed[index].wait(timeout=10)
            return observations

        def run_b():
            observations = []
            for index in range(repetitions):
                assert a_configured[index].wait(timeout=10)
                _configure(dll_b, linelist, case_b)
                b_configured[index].set()
                assert a_observed[index].wait(timeout=10)
                observations.append(_observe(dll_b))
                b_observed[index].set()
            return observations

        with ThreadPoolExecutor(max_workers=2) as executor:
            future_a = executor.submit(run_a)
            future_b = executor.submit(run_b)
            observations_a = future_a.result(timeout=60)
            observations_b = future_b.result(timeout=60)
    finally:
        dll_a.SetEosWarmStartMode(False)

    mismatches = {
        "A": _mismatch_counts(observations_a, reference_a),
        "B": _mismatch_counts(observations_b, reference_b),
    }
    summary = ", ".join(
        f"{instance}.{name}={count}/{repetitions}"
        for instance, instance_counts in mismatches.items()
        for name, count in instance_counts.items()
    )
    assert not any(
        count
        for instance_counts in mismatches.values()
        for count in instance_counts.values()
    ), summary


@skipif_smelib
def test_process_lock_preserves_flux_and_eos_outputs():
    linelist = _make_linelist()
    case_a, case_b = _make_cases()
    dll_a = SME_DLL()
    dll_b = SME_DLL()
    dll_a.SetEosWarmStartMode(True)

    try:
        with dll_a.session():
            reference_a = _run_case(dll_a, linelist, case_a)
        with dll_b.session():
            reference_b = _run_case(dll_b, linelist, case_b)

        repetitions = 20
        start_iteration = Barrier(2, timeout=10)

        def run(dll, case):
            observations = []
            for _ in range(repetitions):
                start_iteration.wait()
                with dll.session():
                    observations.append(_run_case(dll, linelist, case))
            return observations

        with ThreadPoolExecutor(max_workers=2) as executor:
            future_a = executor.submit(run, dll_a, case_a)
            future_b = executor.submit(run, dll_b, case_b)
            observations_a = future_a.result(timeout=60)
            observations_b = future_b.result(timeout=60)
    finally:
        dll_a.SetEosWarmStartMode(False)

    for actual in observations_a:
        _assert_observation_matches(actual, reference_a)
    for actual in observations_b:
        _assert_observation_matches(actual, reference_b)
