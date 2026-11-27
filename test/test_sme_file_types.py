"""Loading an SME file must not import or call anything the file names."""

import io
import json
import os
import sys
import tarfile
import zipfile
from datetime import datetime

import numpy as np
import pandas as pd
import pytest
from astropy import units as u
from astropy.coordinates import (
    ICRS,
    Angle,
    Distance,
    EarthLocation,
    Latitude,
    Longitude,
    SkyCoord,
    frame_transform_graph,
)
from astropy.utils import data as astropy_data
from astropy.utils import iers

# PySME first: it restores two numpy functions that flex-format imports.
from pysme.sme import SME_Structure  # isort: skip
from flex.flex import FlexFile  # isort: skip
from flex.json_encoder import FlexJSONDecoder  # isort: skip
from pysme import sme_reader

load = SME_Structure.load


def structure():
    sme = SME_Structure()
    sme.teff, sme.logg = 5700, 4.4
    sme.wave = np.linspace(6000.0, 6001.0, 50)
    sme.spec = np.linspace(0.9, 1.0, 50)
    return sme


def rewrite(source, target, change=None, add=None):
    """Copy an SME file (a tar); change the JSON of members, or add members."""
    with tarfile.open(source) as old, tarfile.open(target, "w") as new:
        names = old.getnames()
        assert not change or set(change) <= set(names), names
        members = [(i.name, old.extractfile(i).read()) for i in old.getmembers()]
        for name, data in members + list((add or {}).items()):
            if change and name in change:
                tree = json.loads(data)
                change[name](tree)
                data = json.dumps(tree).encode()
            info = tarfile.TarInfo(name)
            info.size = len(data)
            new.addfile(info, io.BytesIO(data))


@pytest.fixture
def saved(tmp_path):
    path = str(tmp_path / "plain.sme")
    structure().save(path)
    return path


def asks_for_a_call(marker):
    # Harmless: all that loading this could do is make an empty folder.
    return {"__module__": "os", "__class__": "mkdir", "value": str(marker)}


def refused(path, saying):
    """Load a file that must be refused; return only after checking why."""
    try:
        load(path)
    except Exception as exc:
        error = exc
    else:
        error = None
    # What matters most is checked by the caller: that nothing was done.
    return lambda: (
        type(error).__name__ == "UnsafeFileError" and saying in str(error)
    ) or pytest.fail(f"not refused as expected: {error!r}")


@pytest.mark.parametrize(
    "member, put",
    [
        ("header.json", lambda tree, value: tree.update(note=value)),
        ("abund/header.json", lambda tree, value: tree.update(note=value)),
        (
            "system_info/header.json",
            lambda tree, value: tree.update(deeper={"a": [{"b": value}]}),
        ),
    ],
)
def test_a_value_naming_a_function_does_not_call_it(tmp_path, saved, member, put):
    marker, changed = tmp_path / "called", str(tmp_path / "changed.sme")
    value = asks_for_a_call(marker)
    rewrite(saved, changed, change={member: lambda tree: put(tree, value)})
    check = refused(changed, "os.mkdir")
    assert not marker.exists()
    check()


def test_a_value_inside_a_plain_json_part_does_not_call_it(tmp_path, saved):
    marker, changed = tmp_path / "called", str(tmp_path / "changed.sme")
    header = {
        "__module__": "flex.extensions.jsondata",
        "__class__": "JsonDataExtension",
        "__header__": True,
    }
    rewrite(
        saved,
        changed,
        add={
            "extra/header.json": json.dumps(header).encode(),
            "extra/data.json": json.dumps({"x": [asks_for_a_call(marker)]}).encode(),
        },
    )
    check = refused(changed, "os.mkdir")
    assert not marker.exists()
    check()


def test_a_part_naming_another_module_does_not_get_it_imported(
    tmp_path, saved, monkeypatch
):
    marker = tmp_path / "imported"
    (tmp_path / "planted_module.py").write_text(
        f"import os\nos.mkdir({str(marker)!r})\nclass Thing:\n    pass\n"
    )
    monkeypatch.syspath_prepend(str(tmp_path))
    changed = str(tmp_path / "changed.sme")
    rename = lambda tree: tree.update(__module__="planted_module", __class__="Thing")
    rewrite(saved, changed, change={"abund/header.json": rename})
    check = refused(changed, "planted_module.Thing")
    assert not marker.exists() and "planted_module" not in sys.modules
    check()


@pytest.mark.parametrize(
    "module, name",
    [(["os"], "mkdir"), ("flex.extensions.tabledata", "TableExtension")],
)
def test_other_names_and_other_flex_parts_are_refused(tmp_path, saved, module, name):
    changed = str(tmp_path / "changed.sme")
    rename = lambda tree: tree.update(__module__=module, __class__=name)
    rewrite(saved, changed, change={"abund/header.json": rename})
    refused(changed, "asks for the Python type")()


def test_flex_formats_own_reader_takes_no_part(saved, monkeypatch):
    def used(*args, **kwargs):
        raise AssertionError("flex-format's reader was used")

    monkeypatch.setattr(FlexFile, "read", used)
    monkeypatch.setattr(FlexFile, "_read_ext_class", used)
    monkeypatch.setattr(FlexJSONDecoder, "_object_hook", used)
    assert load(saved).teff == 5700


def test_what_pysme_and_astropy_store_comes_back(tmp_path):
    sme = structure()
    # (flex-format's writer cannot save a coordinate in a frame without ra and
    # dec, and of a SkyCoord it keeps only the frame, ra and dec, not distance or
    # velocities. astropy Time and SpectralCoord are deliberately not read.)
    sme.meta = {
        "made": datetime(2026, 1, 2, 10, 30),
        "unit": u.erg / u.s / u.cm**2,
        "speeds": u.Quantity([1.0, 2.5], u.km / u.s),
        "angle": Angle(12 * u.deg),
        "longitude": Longitude(370 * u.deg),
        "latitude": Latitude(-5 * u.deg),
        "distance": Distance(10 * u.pc),
        "magnitude": u.Magnitude(12.3),
        "dex": u.Dex(0.5),
        "decibel": u.Decibel(3.0),
        "place": EarthLocation(1 * u.m, 2 * u.m, 3 * u.m),
        "position": SkyCoord(10 * u.deg, 20 * u.deg),
    }
    path = str(tmp_path / "normal.sme")
    sme.save(path)
    loaded = load(path)
    assert (loaded.teff, loaded.logg) == (5700, 4.4)
    assert np.array_equal(np.ravel(loaded.wave), np.ravel(sme.wave))
    for key, value in sme.meta.items():
        assert type(loaded.meta[key]) is type(value), key
        assert repr(loaded.meta[key]) == repr(value), key


@pytest.mark.parametrize(
    "value, change, saying",
    [
        (1.5 * u.km, {"dtype": "int64"}, "does not write"),
        (u.km, {"parse_strict": "silent"}, "does not write"),
        (EarthLocation(1 * u.m, 2 * u.m, 3 * u.m), {"ellipsoid": "WGS72"}, "not write"),
        (SkyCoord(10 * u.deg, 20 * u.deg), {"obstime": "2026-01-01"}, "not write"),
        (SkyCoord(10 * u.deg, 20 * u.deg), {"frame": "altaz"}, "frame"),
        (datetime(2026, 1, 2), {"tzinfo": "UTC"}, "time zone"),
    ],
)
def test_a_stored_value_gets_only_what_pysme_would_write(
    tmp_path, value, change, saying
):
    sme = structure()
    sme.meta = {"x": value}
    path, changed = str(tmp_path / "value.sme"), str(tmp_path / "changed.sme")
    sme.save(path)
    assert repr(load(path).meta["x"]) == repr(value)
    rewrite(path, changed, change={"header.json": lambda t: t["meta"]["x"].update(change)})
    refused(changed, saying)()


def stored_time(**fields):
    return {"__module__": "astropy.time.core", "__class__": "Time", **fields}


def test_a_stored_time_is_refused_and_nothing_is_looked_up(
    tmp_path, saved, monkeypatch
):
    asked = []

    def lookup(*args, **kwargs):
        asked.append(1)
        raise OSError("no network in this test")

    monkeypatch.setattr(iers.IERS_Auto, "open", classmethod(lambda c, *a: lookup()))
    monkeypatch.setattr(astropy_data, "download_file", lookup)
    utc = stored_time(format="isot", value="2026-01-02T00:00:00.000")
    ut1 = stored_time(format="fits", value="2026-01-02T00:00:00.000(UT1)")
    changed = str(tmp_path / "changed.sme")
    # two time scales in one value made astropy fetch the IERS tables
    for value in (stored_time(value=[utc, ut1]), ut1, utc):
        put = lambda tree, value=value: tree.update(meta={"t": value})
        rewrite(saved, changed, change={"header.json": put})
        check = refused(changed, "astropy.time.core.Time")
        assert not asked
        check()


@pytest.fixture
def position(tmp_path):
    sme = structure()
    sme.meta = {"x": SkyCoord(10 * u.deg, 20 * u.deg)}
    path = str(tmp_path / "position.sme")
    sme.save(path)
    return path


def test_a_stored_position_must_name_its_frame(tmp_path, position):
    changed = str(tmp_path / "changed.sme")
    drop = lambda tree: tree["meta"]["x"].pop("frame")
    rewrite(position, changed, change={"header.json": drop})
    refused(changed, "frame")()


def test_the_frame_is_not_looked_up_by_its_name(position, monkeypatch):
    # as if a program had put a frame of its own under an astropy name
    class Impostor(ICRS):
        pass

    lookup_name = frame_transform_graph.lookup_name
    monkeypatch.setattr(
        frame_transform_graph,
        "lookup_name",
        lambda name: Impostor if name == "icrs" else lookup_name(name),
    )
    assert type(load(position).meta["x"].frame) is ICRS


def test_members_are_read_one_at_a_time_and_only_when_needed(
    tmp_path, saved, monkeypatch
):
    changed = str(tmp_path / "changed.sme")
    rewrite(saved, changed, add={"abund/not_needed.bin": b"x" * 100})
    events = []
    extract, numpy_load = tarfile.TarFile.extractfile, np.load

    def extractfile(self, member):
        events.append(getattr(member, "name", member))
        return extract(self, member)

    def decode(*args, **kwargs):
        events.append("decoded")
        return numpy_load(*args, **kwargs)

    monkeypatch.setattr(tarfile.TarFile, "extractfile", extractfile)
    monkeypatch.setattr(np, "load", decode)
    assert load(changed).teff == 5700
    assert "abund/not_needed.bin" not in events
    arrays = [i for i, name in enumerate(events) if name.endswith(".npy")]
    assert len(arrays) > 2
    # each array is decoded before the next member is read
    assert all(events[i + 1] == "decoded" for i in arrays), events


def old_format(path, table_member, table):
    """A file in the zip format of early PySME versions."""
    with zipfile.ZipFile(path, "w") as archive:
        archive.writestr("info.json", json.dumps({"teff": 5500}))
        archive.writestr(
            "linelist/info.json",
            json.dumps({"format": "short", "medium": "air", "citation_info": ""}),
        )
        archive.writestr(f"linelist/{table_member}", table)


LINES = json.dumps(
    [
        {
            "species": "Fe 1",
            "wlcent": 6000.5,
            "gflog": -1.0,
            "excit": 2.0,
            "vmic": 0.0,
            "atom_number": 26.0,
            "ionization": 1.0,
            "gamrad": 8.0,
            "gamqst": -6.0,
            "gamvw": -7.5,
            "lande": 1.0,
            "depth": 0.5,
            "reference": "none",
        }
    ]
)


def test_an_old_format_file_still_opens(tmp_path):
    path = str(tmp_path / "old.sme")
    old_format(path, "data.json", LINES)
    loaded = load(path)
    assert loaded.teff == 5500 and len(loaded.linelist) == 1


def test_an_old_feather_line_list_is_not_given_to_an_old_reader(tmp_path, monkeypatch):
    read = []
    monkeypatch.setattr(pd, "read_feather", lambda *a, **k: read.append(1))
    monkeypatch.setattr(sme_reader, "version", lambda name: "21.0.0")
    path = str(tmp_path / "old.sme")
    old_format(path, "data.feather", b"not really feather")
    check = refused(path, "Feather")
    assert not read
    check()


@pytest.mark.parametrize(
    "installed, safe",
    [
        (None, False),
        ("0.14.0", False),
        ("14.0.0", False),
        ("14.0.1rc1", False),
        ("14.0.1.dev3", False),
        ("14.0.1", True),
        ("14.0.1.post1", True),
        ("21.0.0", True),
    ],
)
def test_trusted_feather_still_needs_a_recent_final_pyarrow(
    monkeypatch, installed, safe
):
    def version(name):
        if installed is None:
            raise sme_reader.PackageNotFoundError(name)
        return installed

    monkeypatch.setattr(sme_reader, "version", version)
    with sme_reader.trusting(True):
        if safe:
            sme_reader.require_safe_feather()
        else:
            with pytest.raises(sme_reader.UnsafeFileError, match="14.0.1"):
                sme_reader.require_safe_feather()
    # whatever pyarrow there is, an untrusted file is refused
    with pytest.raises(sme_reader.UnsafeFileError, match="trust"):
        sme_reader.require_safe_feather()


def test_a_trusted_feather_file_is_read_and_the_trust_ends_with_it(
    tmp_path, monkeypatch
):
    read = []
    table = pd.read_json(io.StringIO(LINES), orient="records")
    monkeypatch.setattr(pd, "read_feather", lambda *a, **k: read.append(1) or table)
    monkeypatch.setattr(sme_reader, "version", lambda name: "21.0.0")
    path = str(tmp_path / "old.sme")
    old_format(path, "data.feather", b"stands in for feather data")
    check = refused(path, "trust")
    assert not read
    check()
    loaded = SME_Structure.load(path, allow_pickle=True)
    assert len(loaded.linelist) == 1 and read == [1]
    check = refused(path, "trust")
    assert read == [1]
    check()


@pytest.mark.skipif(
    sys.platform == "win32", reason="Windows does not let an open file be replaced"
)
def test_trust_stays_with_the_file_that_was_opened(tmp_path, monkeypatch):
    path, other = tmp_path / "old.sme", tmp_path / "other.sme"
    old_format(str(path), "data.json", LINES)
    old_format(str(other), "data.feather", b"not really feather")
    read = []
    monkeypatch.setattr(pd, "read_feather", lambda *a, **k: read.append(1))
    monkeypatch.setattr(sme_reader, "version", lambda name: "21.0.0")
    read_tar = sme_reader.read

    def swapped_after_the_first_try(source):
        try:
            return read_tar(source)  # fails: this is the old zip format
        finally:
            if other.exists():
                os.replace(other, path)

    monkeypatch.setattr(sme_reader, "read", swapped_after_the_first_try)
    loaded = SME_Structure.load(str(path), allow_pickle=True)
    assert not other.exists()  # the swap did happen
    assert not read and len(loaded.linelist) == 1


@pytest.mark.parametrize("suffix", [".npy", ".npz"])
def test_pickle_files_are_loaded_only_when_asked_for(tmp_path, monkeypatch, suffix):
    calls = []

    def numpy_load(path, **kwargs):
        calls.append(kwargs.get("allow_pickle"))
        raise OSError("stop here")

    monkeypatch.setattr(np, "load", numpy_load)
    path = str(tmp_path / ("file" + suffix))
    open(path, "wb").close()
    check = refused(path, "pickle")
    assert calls == []
    check()
    with pytest.raises(OSError, match="stop here"):
        SME_Structure.load(path, allow_pickle=True)
    assert calls == [True]
