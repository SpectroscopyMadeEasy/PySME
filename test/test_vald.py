# -*- coding: utf-8 -*-
from os.path import dirname, join

import numpy as np
import pandas as pd
import pytest

from pysme.abund import Abund
from pysme.linelist.linelist import LineList
from pysme.linelist.vald import ValdFile
from pysme.sme import SME_Structure

species = "Fe 1"
wlcent = 5502.9931
excit = 0.9582
gflog = -3.047
gamrad = 7.19
gamqst = -6.22
gamvw = 239.249
linedata = [species, wlcent, excit, gflog, gamrad, gamqst, gamvw]


def test_line_init():
    """Test that property values equal line data passed to __init__()."""
    line = LineList(linedata)
    assert isinstance(line, LineList)
    assert line.species[0] == species
    assert line.wlcent[0] == wlcent
    assert line.excit[0] == excit
    assert line.gflog[0] == gflog
    assert line.gamrad[0] == gamrad
    assert line.gamqst[0] == gamqst
    assert line.gamvw[0] == gamvw


def test_linelist_add_and_len():
    """Test that len() returns the number of lines (including 0) in list."""
    linelist = LineList()
    assert isinstance(linelist, LineList)
    assert len(linelist) == 0
    for iline in range(3):
        print(len(linelist._lines))
        print(linelist._lines)
        assert len(linelist) == iline
        linelist.add(*linedata)


def test_linelist_properties():
    """Test that properties are lists with one item per line.
    Test that property value equal line data passed to add().
    """
    linelist = LineList()
    linelist.add(*linedata)
    proplist = [
        linelist.species,
        linelist.wlcent,
        linelist.excit,
        linelist.gflog,
        linelist.gamrad,
        linelist.gamqst,
        linelist.gamvw,
    ]
    for iprop, prop in enumerate(proplist):
        print(linelist._lines)
        assert isinstance(prop, np.ndarray)
        assert len(prop) == 1
        assert prop[0] == linedata[iprop]


def test_valdfile():
    """Test class to read a VALD line data file."""
    testdir = dirname(__file__)
    linelist = ValdFile(join(testdir, "testcase1.lin"))

    assert len(linelist) == 44
    assert linelist.lineformat == "short"
    assert linelist.medium == "air"
    assert linelist[0].species == "V 1"

    with pytest.raises(IOError):
        ValdFile(join(testdir, "testcase1.npy"))


def test_medium():
    """Test class to read a VALD line data file."""
    testdir = dirname(__file__)
    vf = ValdFile(join(testdir, "testcase1.lin"), medium="vac")

    assert vf.medium == "vac"


def test_header_without_vmicro_separator():
    """Accept VALD-compatible headers missing the comma after Vmicro."""
    vald = ValdFile.__new__(ValdFile)
    nlines = vald.parse_header(
        "3600.00000, 8000.00000, 25488, 303039, 2.0 Wavelength region, "
        "lines selected, lines processed, Vmicro"
    )

    assert nlines == 25488
    assert vald._nlines_proc == 303039
    assert vald._vmicro == 2.0


def test_short_format():
    linelist = ValdFile(join(dirname(__file__), "testcase1.lin"))

    assert linelist.lineformat == "short"
    assert len(linelist) == 44
    assert linelist.atomic is not None
    assert linelist.species is not None

    with pytest.raises(AttributeError):
        _ = linelist.lulande

    with pytest.raises(AttributeError):
        _ = linelist.extra

    assert isinstance(linelist.abund, Abund)
    assert isinstance(linelist.atmo, str)


def test_long_format():
    linelist = ValdFile(join(dirname(__file__), "testcase3.lin"))

    assert linelist.lineformat == "long"
    assert len(linelist) == 67
    assert linelist.atomic is not None
    assert linelist.species is not None
    assert linelist.lulande is not None
    assert linelist.extra is not None

    assert linelist.reference is not None
    assert linelist.error is not None
    assert linelist.term_lower is not None
    assert linelist.term_upper is not None

    assert isinstance(linelist.abund, Abund)
    assert isinstance(linelist.atmo, str)


def test_streaming_long_parser_matches_legacy_across_chunks():
    path = join(dirname(__file__), "testcase3.lin")
    legacy = ValdFile.__new__(ValdFile)
    legacy.citation_info = ValdFile.citation_info
    legacy_lines = legacy._loads_legacy(path)
    streaming = ValdFile.__new__(ValdFile)
    streaming.citation_info = ValdFile.citation_info
    streaming_lines = streaming._loads_long_extract_stellar_stream(
        path, chunk_size=3
    )

    pd.testing.assert_frame_equal(streaming_lines, legacy_lines)
    assert streaming.header == legacy.header
    assert streaming.valdtype == legacy.valdtype
    assert streaming.lineformat == legacy.lineformat
    assert streaming.medium == legacy.medium
    assert streaming.unit == legacy.unit
    assert streaming.energy_unit == legacy.energy_unit
    assert streaming.atmo == legacy.atmo
    assert streaming._ref_record == legacy._ref_record
    assert streaming.citation_info == legacy.citation_info
    assert streaming.abund.monh == legacy.abund.monh
    assert streaming.abund.type == legacy.abund.type
    assert np.array_equal(streaming.abund._pattern, legacy.abund._pattern)


@pytest.mark.parametrize("fmt, idiscard", [("long", 45), ("short", 90)])
def test_extract_hyphenated_reference_keys(fmt, idiscard):
    # The leading part of a line (data or description) is never searched
    head = "x 1 NOT-A-KEY ".ljust(idiscard)
    lines = [
        head + "   9 wl:K13  10 gf:RU   6 BA-J   3 Si2-av1   3 gf:LD-HS'\n",
        # Several keys joined by "+", and a hyphenated key before the quote
        head + "   2 wl:S-G+K13   4 B-WPNP'\n",
        # Dashes in free text must not produce keys
        head + "   5 K14   no-hfs Fe 7 -x 7 y- 7 z--w '\n",
    ]

    keys = ValdFile._extract_reference_keys(lines, fmt)

    assert keys == {
        "K13",
        "RU",
        "BA-J",
        "Si2-av1",
        "LD-HS",
        "S-G+K13",
        "B-WPNP",
        "K14",
    }


def test_citation_includes_hyphenated_key_long_format():
    linelist = ValdFile(join(dirname(__file__), "testcase3.lin"))

    assert "@article{BA-J," in linelist.citation_info
    sme = SME_Structure()
    sme.linelist = linelist
    assert "Aspelund-Johansson" in sme.citation()


def test_citation_includes_hyphenated_key_short_format(tmp_path):
    with open(join(dirname(__file__), "testcase1.lin")) as file:
        content = file.read()
    assert "4 CB" in content
    path = tmp_path / "hyphen_short.lin"
    path.write_text(content.replace("4 CB", "4 BA-J", 1))

    linelist = ValdFile(path)

    assert linelist.lineformat == "short"
    assert "@article{BA-J," in linelist.citation_info
    sme = SME_Structure()
    sme.linelist = linelist
    assert "Aspelund-Johansson" in sme.citation()


@pytest.mark.parametrize(
    "line, expected",
    [
        # Energies in eV: the entries start well before column 45
        (
            "'_          Kurucz Fe I 2014   1 wl:AA   2 gf:BB   3 CC   4 DD    Fe'",
            {"AA", "BB", "CC", "DD"},
        ),
        # Several spaces between number, tag and key
        (
            "'N C        Kurucz Fe I 2014   1  wl:AA   2  gf:BB   3     CC   no-hfs Fe'",
            {"AA", "BB", "CC"},
        ),
        # A single space between entries
        (
            "'NA         CNO, Na 1, Mg 1:   00 wl:AA 00 gf:BB 00 CC 665 DD 00 EE hfs:XX,YY Na'",
            {"AA", "BB", "CC", "DD", "EE"},
        ),
        # Energies in cm^-1: wider description, with digits in it
        (
            "'N D        Kurucz Fe I 2014Fe               1 wl:AA   2 gf:BB   1 CC'",
            {"AA", "BB", "CC"},
        ),
        # Several keys in one entry, with or without their own tag
        (
            "'_          Sc-Ni Wisconsin   53 wl:AA,iso:BB 41 iso:CC,DD  21 gf:EE,FF-G   4 HH    Ni'",
            {"AA", "BB", "CC", "DD", "EE", "FF-G", "HH"},
        ),
        # Numbers followed by words in the description and the comment
        (
            "'_          La1 to Lu 2 Bell  19 wl:AA  19 gf:BB  19 CC    hfs:XX,YY Pr+'",
            {"AA", "BB", "CC"},
        ),
        (
            "'_          Kurucz MnI 2007 Mn    8 wl:AA   8 gf:BB   8 CC \\\\hfs XX XX'",
            {"AA", "BB", "CC"},
        ),
        # Untagged wavelength entry in front of "gf:"
        (
            "'_          La1 to Lu 2 Bell  19 AA  19 gf:BB,CC  19 DD    La'",
            {"AA", "BB", "CC", "DD"},
        ),
        (
            "'_          Kurucz Fe I 2014   1 gf:AA   1 BB    Fe'",
            {"AA", "BB"},
        ),
    ],
)
def test_extract_reference_keys_long_layouts(line, expected):
    # Repeat the line to check that nothing is read across line breaks
    keys = ValdFile._extract_reference_keys([line + "\n"] * 2, "long")

    assert keys == expected


@pytest.mark.parametrize("padding", [0, 40])
def test_extract_reference_keys_short_layouts(padding):
    data = "'Fe 1', 5000.0, 1.0, 1.0, -1.0, 8.0,-6.0, -7.0, 1.0, 0.5, ".ljust(
        89 + padding - 4
    )
    line = data + "'   1 wl:AA   2 BB   3 gf:CC   4 DD Fe      '\n"

    keys = ValdFile._extract_reference_keys([line], "short")

    assert keys == {"AA", "BB", "CC", "DD"}


def test_extract_reference_keys_without_tags():
    # Without a tagged entry only the part behind the fixed column is searched
    line = "'       K07            Kurucz FeI  Fe".ljust(50) + "1 AA   1 BB '\n"

    keys = ValdFile._extract_reference_keys([line] * 2, "long")

    assert keys == {"AA", "BB"}


@pytest.mark.parametrize(
    "filename, expected",
    [
        ("extract_element.lin", {"BPM", "FMW", "K14"}),
        # "BSSc" and "LGWS" are cut to four characters in the file itself
        (
            "extract_stellar.lin",
            {"BA-J", "BPM", "BSSc", "FMW", "K08", "K13", "K14", "K16"}
            | {"LBS", "LGWS", "NWL", "TB"},
        ),
    ],
)
def test_reference_keys_of_compact_files(filename, expected):
    with open(join(dirname(__file__), filename)) as file:
        lines = file.readlines()
    # Both files are in long format with four text lines per spectral line
    references = [line for line in lines if " wl:" in line]

    assert ValdFile._extract_reference_keys(references, "long") == expected
    citation_info = ValdFile(join(dirname(__file__), filename)).citation_info
    assert "@article{FMW," in citation_info


def test_save_reload_long_subset_updates_header(tmp_path):
    linelist = ValdFile(join(dirname(__file__), "testcase3.lin"))
    subset = linelist[:5]

    path = tmp_path / "subset_long.lin"
    subset.save(path, overwrite=True)
    assert path.read_text().splitlines()[0].split(",")[2].strip() == "5"
    assert "species_len" not in subset._lines.columns
    reloaded = ValdFile(path)

    assert len(reloaded) == len(subset)
    assert reloaded.lineformat == "long"
    assert np.allclose(reloaded.wlcent, subset.wlcent)


def test_save_reload_short_subset_updates_header(tmp_path):
    linelist = ValdFile(join(dirname(__file__), "testcase1.lin"))
    subset = linelist[:5]

    path = tmp_path / "subset_short.lin"
    subset.save(path, overwrite=True)
    assert path.read_text().splitlines()[0].split(",")[2].strip() == "5"
    assert "species_len" not in subset._lines.columns
    reloaded = ValdFile(path)

    assert len(reloaded) == len(subset)
    assert reloaded.lineformat == "short"
    assert np.allclose(reloaded.wlcent, subset.wlcent)
