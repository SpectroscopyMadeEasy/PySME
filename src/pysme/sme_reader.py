# -*- coding: utf-8 -*-
"""
Read SME files as data.

An SME file is a tar archive written with flex-format. flex-format's own reader
imports the modules and calls the classes that a file names, so a file from
someone else could run code. This reader does not use it: it reads the members
itself and builds only a fixed set of types, from the tables below. Nothing is
imported or looked up by a name that comes from the file.

Some values that flex-format could save are therefore not read back, and a file
holding one is refused: among them astropy Time (building one from a file can
make astropy fetch the IERS tables) and SpectralCoord. PySME itself stores
neither; they can only come from a user's own entries, for example in ``meta``.

PySME-GUI carries a copy of this file, for use with older PySME versions;
keep the two the same.
"""

import contextvars
import datetime
import io
import json
import logging
import re
import tarfile
from contextlib import contextmanager
from importlib.metadata import PackageNotFoundError, version

import numpy as np
import pandas as pd
from astropy import coordinates, units
from astropy.coordinates import builtin_frames
from astropy.units.function import logarithmic

# First, because PySME restores two numpy functions that flex-format imports.
import pysme  # noqa: F401  # isort: skip
from flex.extensions.bindata import BinaryDataExtension, MultipleDataExtension
from flex.extensions.jsondata import JsonDataExtension
from flex.extensions.tabledata import JSONTableExtension
from flex.flex import FlexFile

logger = logging.getLogger(__name__)

# Reading Feather data from an unknown source is safe only from this pyarrow on.
FEATHER_NEEDS_PYARROW = (14, 0, 1)


class UnsafeFileError(ValueError):
    """The file asks for something that cannot be read as plain data."""


def _refuse(reason):
    raise UnsafeFileError(
        f"This file was not loaded: {reason}. If you saved this yourself, remove it "
        "from the file's information and save again. If you did not, be careful: "
        "opened with other tools, a file made for it can run programs on your "
        "computer."
    )


def _refuse_type(module, name):
    asked = f"{module}.{name}"[:100]
    _refuse(f"it asks for the Python type '{asked}', which is not read as data")


def _lookup(table, module, name):
    found = None
    if isinstance(module, str) and isinstance(name, str):
        found = table.get((module, name))
    if found is None:
        _refuse_type(module, name)
    return found


# Values that flex-format stores together with their type: dates and
# astronomical quantities. For each, the fields flex-format's writer puts in the
# file; no other field is handed on to the type.
_QUANTITY = ("value", "unit")
_DATETIME = (
    "tzinfo",
    "year",
    "month",
    "day",
    "hour",
    "minute",
    "second",
    "microsecond",
)
_VALUES = {
    ("datetime", "datetime"): (datetime.datetime, _DATETIME),
    ("astropy.units", "Unit"): (units.Unit, ("value",)),
    ("astropy.units.quantity", "Quantity"): (units.Quantity, _QUANTITY),
    ("astropy.coordinates.sky_coordinate", "SkyCoord"): (
        coordinates.SkyCoord,
        ("frame", "ra", "dec"),
    ),
    ("astropy.coordinates.earth", "EarthLocation"): (
        coordinates.EarthLocation,
        ("x", "y", "z"),
    ),
    ("astropy.coordinates.distances", "Distance"): (coordinates.Distance, _QUANTITY),
}
for _cls in (logarithmic.Magnitude, logarithmic.Dex, logarithmic.Decibel):
    _VALUES[(_cls.__module__, _cls.__name__)] = (_cls, _QUANTITY)
# astropy has had the angles in two places.
for _module in ("astropy.coordinates.angles", "astropy.coordinates.angles.core"):
    for _cls in (coordinates.Angle, coordinates.Longitude, coordinates.Latitude):
        _VALUES[(_module, _cls.__name__)] = (_cls, _QUANTITY)

# astropy's own frames that have ra and dec. A file names one; the class comes
# from this table, never from astropy's registry, which a program can add to.
_FRAMES = {
    frame.name: frame
    for frame in (
        builtin_frames.ICRS,
        builtin_frames.FK4,
        builtin_frames.FK4NoETerms,
        builtin_frames.FK5,
        builtin_frames.GCRS,
        builtin_frames.CIRS,
        builtin_frames.HCRS,
        builtin_frames.LSR,
        builtin_frames.LSRD,
        builtin_frames.LSRK,
        builtin_frames.PrecessedGeocentric,
        builtin_frames.TETE,
    )
}


def _check_fields(cls, fields, allowed):
    extra = sorted(set(fields) - set(allowed))
    if extra:
        _refuse(
            f"a stored {cls.__name__} has fields PySME does not write "
            f"({', '.join(extra)[:80]})"
        )
    if cls is coordinates.SkyCoord:
        # Always written; without it astropy would choose one.
        frame = fields.get("frame")
        if not (isinstance(frame, str) and frame in _FRAMES):
            _refuse(f"a stored SkyCoord asks for the frame '{str(frame)[:40]}'")
    if "unit" in fields and not isinstance(fields["unit"], str):
        _refuse(f"the unit of a stored {cls.__name__} is not text")
    if fields.get("tzinfo") is not None:
        _refuse("a stored date asks for a time zone")


def _value(obj):
    """Build a value of a listed type from a JSON object that names one."""
    # The condition under which flex-format builds something from an object.
    if (
        "__class__" not in obj
        or "__module__" not in obj
        or obj.get("__header__", False)
    ):
        return obj
    cls, allowed = _lookup(_VALUES, obj["__module__"], obj["__class__"])
    fields = {k: v for k, v in obj.items() if k not in ("__class__", "__module__")}
    _check_fields(cls, fields, allowed)
    if cls is coordinates.SkyCoord:
        fields["frame"] = _FRAMES[fields["frame"]]
    if cls is units.Quantity and fields.get("value") is None:
        fields["value"] = float("nan")
    args = (fields.pop("value"),) if "value" in fields else ()
    try:
        return cls(*args, **fields)
    except Exception as ex:
        logger.warning("Could not rebuild a %s: %s", cls.__name__, ex)
        return obj


def _json(data):
    return json.loads(data, object_hook=_value)


def _array(data):
    return np.load(io.BytesIO(data), allow_pickle=False)


# The parts a file is made of. PySME writes the two binary ones and the JSON
# table; the plain JSON part is flex-format's own. Each builder is given the
# names of the part's members and a function that reads one of them, so that only
# what it needs is read, one member at a time.
def _binary(header, members, content):
    return BinaryDataExtension(header, _array(content(members["data.npy"])))


def _multiple(header, members, content):
    data = {name[:-4]: _array(content(full)) for name, full in members.items()}
    return MultipleDataExtension(header, data=data)


def _plain_json(header, members, content):
    data = _json(content(members["data.json"])) if "data.json" in members else {}
    return JsonDataExtension(header, data)


def _json_table(header, members, content):
    text = io.BytesIO(content(members["data.json"]))
    return JSONTableExtension(header, pd.read_json(text, orient="records"))


_PARTS = {
    ("flex.extensions.bindata", "BinaryDataExtension"): _binary,
    ("flex.extensions.bindata", "MultipleDataExtension"): _multiple,
    ("flex.extensions.jsondata", "JsonDataExtension"): _plain_json,
    ("flex.extensions.tabledata", "JSONTableExtension"): _json_table,
}


def read(filename):
    """
    Read an SME file in the current (tar) format

    Parameters
    ----------
    filename : str or binary file
        file to read; an open file is left open

    Returns
    -------
    FlexFile
        the same header and parts that flex-format would give

    Raises
    ------
    UnsafeFileError
        if the file names a type that is not on the lists above
    """
    if hasattr(filename, "read"):
        opened = tarfile.open(fileobj=filename, mode="r")
    else:
        opened = tarfile.open(filename, mode="r")
    with opened as archive:

        def content(name):
            member = archive.extractfile(name)
            if member is None:
                raise ValueError(f"{name!r} in the SME file is not a file")
            return member.read()

        header = _json(content("header.json"))
        # The first part of a member's name says which part it belongs to.
        groups = {}
        for name in archive.getnames():
            if name != "header.json":
                group = name.replace("\\", "/").partition("/")[0]
                groups.setdefault(group, []).append(name)

        extensions = {}
        for group in sorted(groups):
            headers = [n for n in groups[group] if n.endswith("header.json")]
            if not headers:
                raise ValueError(f"The part {group!r} of the SME file has no header")
            part = _json(content(headers[-1]))
            fields = part if isinstance(part, dict) else {}
            build = _lookup(_PARTS, fields.get("__module__"), fields.get("__class__"))
            members = {n[len(group) + 1 :]: n for n in groups[group] if n not in headers}
            extensions[group] = build(part, members, content)

    return FlexFile(header=header, extensions=extensions)


# Set while loading a file the caller has said they trust.
_trusted = contextvars.ContextVar("sme_file_trusted", default=False)


@contextmanager
def trusting(trusted):
    """Within this block, the file being loaded is trusted if ``trusted``."""
    token = _trusted.set(bool(trusted))
    try:
        yield
    finally:
        _trusted.reset(token)


def require_safe_feather():
    """
    Raise unless Feather data may be read: the caller trusts the file, and
    pyarrow is recent enough.

    A recent pyarrow alone is not enough for a file from an unknown source: a
    program can switch the unsafe behaviour back on, or have added decoders of
    its own, which a file can then choose.
    """
    needed = ".".join(map(str, FEATHER_NEEDS_PYARROW))
    if not _trusted.get():
        raise UnsafeFileError(
            "This file was not loaded: it keeps its line list in the Feather format, "
            "as some old SME files do, and reading that format can run code hidden "
            "in a file. If you made the file yourself, or trust where it comes "
            "from, load it with SME_Structure.load(filename, allow_pickle=True); "
            f"this needs pyarrow {needed} or newer."
        )
    try:
        found = version("pyarrow")
    except PackageNotFoundError:
        found = None
    # Final releases only: "14.0.1rc1" comes before 14.0.1.
    numbers = re.fullmatch(r"(\d+)\.(\d+)\.(\d+)(\.post\d+)?(\+.*)?", found or "")
    if not numbers or tuple(map(int, numbers.groups()[:3])) < FEATHER_NEEDS_PYARROW:
        raise UnsafeFileError(
            "This file was not loaded: reading its Feather line list needs pyarrow "
            f"{needed} or newer (found: {found or 'not installed'})."
        )
