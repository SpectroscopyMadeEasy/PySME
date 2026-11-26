# -*- coding: utf-8 -*-
"""Checks for the bundled VALD bibliography, VALD3_ref.bib."""
import re
from os.path import dirname, join

import pybtex
import pybtex.database

from pysme.data_structure import Collection
from pysme.linelist import vald


def load_bibliography():
    return pybtex.database.parse_file(join(dirname(vald.__file__), "VALD3_ref.bib"))


def test_every_reference_renders():
    """One entry that cannot be rendered makes all of sme.citation() fail."""
    bibliography = load_bibliography()
    failures = {}
    for key, entry in bibliography.entries.items():
        citation_info = pybtex.database.BibliographyData({key: entry}).to_string(
            "bibtex"
        )
        try:
            Collection().create_citation(citation_info)
        except pybtex.exceptions.PybtexError as ex:
            failures[key] = str(ex)

    assert len(bibliography.entries) > 0
    assert failures == {}


def test_all_references_render_together():
    citation_info = load_bibliography().to_string("bibtex")

    assert Collection().create_citation(citation_info)


def test_no_field_is_swallowed_by_another():
    """Unbalanced braces can put the following fields into one value."""
    pattern = re.compile(
        r"\b(title|journal|year|doi|adsurl)\s*=\s*[{\"]|VALD_ref_number"
    )
    damaged = []
    for key, entry in load_bibliography().entries.items():
        values = list(entry.fields.values())
        values += [str(p) for persons in entry.persons.values() for p in persons]
        if any(pattern.search(value) for value in values):
            damaged.append(key)

    assert damaged == []
