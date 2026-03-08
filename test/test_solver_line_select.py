# -*- coding: utf-8 -*-
from os.path import dirname, join

import numpy as np
import pytest

from pysme.sme import SME_Structure as SME_Struct
from pysme.solve import solve
from pysme.synthesize import synthesize_spectrum


cwd = dirname(__file__)
filename = join(cwd, "testcase1.inp")


def _concat_segments(vec, nseg):
    return np.concatenate([np.asarray(vec[i], dtype=float) for i in range(nseg)])


def _prepare_sme():
    sme = SME_Struct.load(filename)
    wran = np.asarray(sme.wran, dtype=float).reshape(-1, 2)
    wl = np.asarray(sme.linelist["wlcent"], dtype=float)
    lo = float(np.min(wran[:, 0]) - 1.0)
    hi = float(np.max(wran[:, 1]) + 1.0)
    keep = (wl >= lo) & (wl <= hi)
    if np.any(keep):
        sme.linelist = sme.linelist[keep]
    return sme


def _run_solve(method):
    sme = _prepare_sme()
    sme.line_select_method = method
    sme.line_select_policy = "strict" if method in ("almax", "cdr") else "auto"
    sme.line_select_recompute = "always"
    sme.line_select_reuse = "none"
    sme.line_select_parallel = False

    # Use only vrad to keep runtime low while still exercising the solve path.
    out = solve(sme, ["vrad"], linelist_mode="all")
    wave = _concat_segments(out.wave, out.nseg)
    flux = _concat_segments(out.synth, out.nseg)
    return out, wave, flux


@pytest.mark.parametrize(
    "method,required_cols",
    [
        ("internal", []),
        ("almax", ["almax_ratio", "line_range_s", "line_range_e", "strong"]),
        ("cdr", ["central_depth", "line_range_s", "line_range_e", "strong"]),
    ],
)
def test_solve_line_select_methods_run(method, required_cols):
    out, _, flux = _run_solve(method)
    cols = set(out.linelist._lines.columns)

    assert out.synth is not None
    assert np.all(np.isfinite(flux))
    for col in required_cols:
        assert col in cols


def test_solve_line_select_methods_match_internal_flux():
    _, wave_ref, flux_ref = _run_solve("internal")

    for method in ("almax", "cdr"):
        _, wave, flux = _run_solve(method)
        flux_interp = np.interp(wave_ref, wave, flux)
        diff = flux_interp - flux_ref
        assert np.max(np.abs(diff)) < 1e-10


def test_line_precompute_database_separates_method_and_linelist_hash(tmp_path):
    db = str(tmp_path)

    sme_cdr = _prepare_sme()
    sme_cdr.line_select_method = "cdr"
    sme_cdr.line_select_policy = "strict"
    sme_cdr.line_select_recompute = "always"
    sme_cdr.line_select_parallel = False
    _ = synthesize_spectrum(
        sme_cdr,
        linelist_mode="all",
        line_precompute_database=db,
    )

    cdr_files = list(tmp_path.glob("cdr_*.npz"))
    assert len(cdr_files) > 0

    sme_almax = _prepare_sme()
    sme_almax.line_select_method = "almax"
    sme_almax.line_select_policy = "strict"
    sme_almax.line_select_recompute = "always"
    sme_almax.line_select_parallel = False
    sme_almax.line_select_almax_threshold = sme_almax.accrt
    _ = synthesize_spectrum(
        sme_almax,
        linelist_mode="all",
        line_precompute_database=db,
    )

    almax_files = list(tmp_path.glob("almax_*.npz"))
    assert len(almax_files) > 0

    # Mismatched linelist (different order/content) should be treated as cache miss.
    sme_mismatch = _prepare_sme()
    sme_mismatch.linelist = sme_mismatch.linelist[1:]
    sme_mismatch.line_select_method = "cdr"
    sme_mismatch.line_select_policy = "strict"
    sme_mismatch.line_select_recompute = "never"
    sme_mismatch.line_select_parallel = False
    with pytest.raises(ValueError, match="no matching entry"):
        _ = synthesize_spectrum(
            sme_mismatch,
            linelist_mode="all",
            line_precompute_database=db,
        )
