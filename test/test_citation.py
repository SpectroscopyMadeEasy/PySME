# -*- coding: utf-8 -*-
import re

import numpy as np
import pytest

from pysme.sme import SME_Structure

# The SME and PySME papers, newest first
CODE_DOIS = [
    "10.1051/0004-6361/202659986",  # Jian et al. 2026
    "10.1051/0004-6361/202244482",  # Wehrhahn et al. 2023
    "10.1051/0004-6361/201629124",  # Piskunov & Valenti 2017
]
CODE_KEYS = [
    "2026A&A...711A.168J",
    "2023A&A...671A.171W",
    "2017A&A...597A..16P",
    "1996A&AS..118..595V",
]
CODE_TITLES = [
    "Improved modelling of stellar spectra",
    "Spectroscopy made easier",
    "Spectroscopy made easy: evolution",
    "A new tool for fitting observations",
]

# Default ``citation_info`` stored in files written before the PySME papers
# were added. Shortened fields would no longer match, so only the whitespace
# may differ from what those files hold.
OLD_DEFAULT = r"""
    @ARTICLE{2017A&A...597A..16P,
        author = {{Piskunov}, Nikolai and {Valenti}, Jeff A.},
        title = "{Spectroscopy Made Easy: Evolution}",
        journal = {\aap},
        keywords = {stars: abundances, radiative transfer, stars: fundamental parameters, stars: atmospheres, techniques: spectroscopic, Astrophysics - Instrumentation and Methods for Astrophysics, Astrophysics - Solar and Stellar Astrophysics},
        year = "2017",
        month = "Jan",
        volume = {597},
        eid = {A16},
        pages = {A16},
        doi = {10.1051/0004-6361/201629124},
        archivePrefix = {arXiv},
        eprint = {1606.06073},
        primaryClass = {astro-ph.IM},
        adsurl = {https://ui.adsabs.harvard.edu/abs/2017A&A...597A..16P},
        adsnote = {Provided by the SAO/NASA Astrophysics Data System}
    }
    @ARTICLE{1996A&AS..118..595V,
        author = {{Valenti}, J.~A. and {Piskunov}, N.},
        title = "{Spectroscopy made easy: A new tool for fitting observations with synthetic spectra.}",
        journal = {\aaps},
        keywords = {RADIATIVE TRANSFER, METHODS: NUMERICAL, TECHNIQUES: SPECTROSCOPIC, STARS: FUNDAMENTAL PARAMETERS, SUN: FUNDAMENTAL PARAMETERS, ATOMIC DATA},
        year = "1996",
        month = "Sep",
        volume = {118},
        pages = {595-603},
        adsurl = {https://ui.adsabs.harvard.edu/abs/1996A&AS..118..595V},
        adsnote = {Provided by the SAO/NASA Astrophysics Data System}
    }
"""


def entry_keys(text):
    return re.findall(r"@\w+\{([^,\s]+),", text)


def assert_in_order(text, fragments):
    # Line wrapping and capitalisation depend on the output backend
    text = " ".join(text.split()).casefold()
    positions = [text.find(f.casefold()) for f in fragments]
    assert -1 not in positions, dict(zip(fragments, positions))
    assert positions == sorted(positions), dict(zip(fragments, positions))


def test_bibtex_citation_lists_code_papers_newest_first():
    text = SME_Structure().citation(output="bibtex")

    assert entry_keys(text)[: len(CODE_KEYS)] == CODE_KEYS
    assert_in_order(text, CODE_DOIS)
    for key in CODE_KEYS:
        assert entry_keys(text).count(key) == 1


def test_string_citation_lists_code_papers_newest_first():
    text = SME_Structure().citation()

    assert_in_order(text, CODE_TITLES)
    assert_in_order(text, CODE_DOIS)
    # The code papers come ahead of the data references
    assert text.startswith("[1] ")
    assert_in_order(text, [CODE_TITLES[-1], "The Solar Chemical Composition"])


@pytest.mark.parametrize("extension", ["sme", "npy"])
def test_old_default_citation_is_updated_on_load(tmp_path, extension):
    sme = SME_Structure()
    sme.citation_info = OLD_DEFAULT
    path = str(tmp_path / f"old.{extension}")
    if extension == "sme":
        sme.save(path)
    else:
        np.save(path, np.array([sme], dtype=object))

    # allow_pickle: this test made the file itself (it matters only for .npy).
    loaded = SME_Structure.load(path, allow_pickle=True)

    assert loaded.citation_info == SME_Structure().citation_info
    assert entry_keys(loaded.citation(output="bibtex"))[:4] == CODE_KEYS
    assert_in_order(loaded.citation(), CODE_TITLES)


def test_customised_citation_is_kept_on_load(tmp_path):
    custom = OLD_DEFAULT + "@misc{mine, title={My reference}, year=2001}\n"
    sme = SME_Structure()
    sme.citation_info = custom
    path = str(tmp_path / "custom.sme")
    sme.save(path)

    assert SME_Structure.load(path).citation_info == custom
