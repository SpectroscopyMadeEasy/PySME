# -*- coding: utf-8 -*-
import hashlib
import os
import shutil
import tarfile
from pathlib import Path

import pytest

from pysme import large_file_storage as lfs_module
from pysme.large_file_storage import LargeFileStorage


def test_lfs_storage_tilde_path_is_expanded(monkeypatch, tmp_path):
    fake_home = tmp_path / "home"
    fake_home.mkdir()
    monkeypatch.setenv("HOME", str(fake_home))

    lfs = LargeFileStorage(
        server="http://example.com",
        pointers={},
        storage="~/.sme/test_cache",
    )

    expected = Path(str(fake_home)).resolve() / ".sme" / "test_cache"
    assert lfs.current == expected
    assert os.environ["XDG_CACHE_HOME"] == str(expected)
    assert expected.exists()


def test_get_urls_supports_multiple_servers(tmp_path):
    lfs = LargeFileStorage(
        server=["https://mirror-a.example", "https://mirror-b.example/"],
        pointers={"foo.grd": "nlte_grids/foo_v1.grd.gz"},
        storage=str(tmp_path / "cache"),
    )

    assert lfs.get_urls("foo.grd") == [
        "https://mirror-a.example/nlte_grids/foo_v1.grd.gz",
        "https://mirror-b.example/nlte_grids/foo_v1.grd.gz",
    ]


def test_get_urls_supports_pointer_lists_and_absolute_urls(tmp_path):
    lfs = LargeFileStorage(
        server=["https://mirror-a.example"],
        pointers={
            "foo.grd": [
                "https://zenodo.org/records/1/files/foo.grd?download=1",
                "nlte_grids/foo_v1.grd.gz",
            ]
        },
        storage=str(tmp_path / "cache"),
    )

    assert lfs.get_urls("foo.grd") == [
        "https://zenodo.org/records/1/files/foo.grd?download=1",
        "https://mirror-a.example/nlte_grids/foo_v1.grd.gz",
    ]


def test_get_urls_supports_pointer_dict_metadata(tmp_path):
    lfs = LargeFileStorage(
        server=["https://mirror-a.example"],
        pointers={
            "foo.grd": [
                {
                    "url": "https://zenodo.org/records/1/files/foo.grd?download=1",
                    "md5": "00112233445566778899aabbccddeeff",
                    "sha256": "abc123",
                    "size": 42,
                },
                {
                    "url": "nlte_grids/foo_v1.grd.gz",
                    "md5": "ffeeddccbbaa99887766554433221100",
                    "sha256": "def456",
                },
            ]
        },
        storage=str(tmp_path / "cache"),
    )

    targets = lfs.get_targets("foo.grd")

    assert [target.url for target in targets] == [
        "https://zenodo.org/records/1/files/foo.grd?download=1",
        "https://mirror-a.example/nlte_grids/foo_v1.grd.gz",
    ]
    assert targets[0].md5 == "00112233445566778899aabbccddeeff"
    assert targets[0].sha256 == "abc123"
    assert targets[0].size == 42
    assert targets[1].md5 == "ffeeddccbbaa99887766554433221100"
    assert targets[1].sha256 == "def456"
    assert targets[1].size is None


def test_get_falls_back_to_next_mirror(monkeypatch, tmp_path):
    payload = tmp_path / "payload.bin"
    payload.write_bytes(b"test")

    calls = []

    def fake_download_file(url, cache=True, pkgname=""):
        calls.append(url)
        if "mirror-a" in url:
            raise OSError("mirror-a unreachable")
        return str(payload)

    monkeypatch.setattr(lfs_module, "download_file", fake_download_file)

    lfs = LargeFileStorage(
        server=["https://mirror-a.example", "https://mirror-b.example"],
        pointers={"foo.grd": "nlte_grids/foo_v1.grd"},
        storage=str(tmp_path / "cache"),
    )

    got = lfs.get("foo.grd")

    assert got == str(payload)
    assert calls == [
        "https://mirror-a.example/nlte_grids/foo_v1.grd",
        "https://mirror-b.example/nlte_grids/foo_v1.grd",
    ]


def test_get_falls_back_when_first_payload_is_html_error_page(monkeypatch, tmp_path):
    bad_payload = tmp_path / "bad.html"
    bad_payload.write_text("<!DOCTYPE html><html><body>login required</body></html>")
    good_payload = tmp_path / "good.grd"
    good_payload.write_bytes(b"grid-data")

    calls = []

    def fake_download_file(url, cache=True, pkgname=""):
        calls.append(url)
        if "nadc" in url:
            return str(bad_payload)
        return str(good_payload)

    monkeypatch.setattr(lfs_module, "download_file", fake_download_file)

    lfs = LargeFileStorage(
        server=["https://mirror.example"],
        pointers={
            "foo.grd": [
                "https://nadc.example/res/file_upload/download?id=1",
                "nlte_grids/foo_v1.grd",
            ]
        },
        storage=str(tmp_path / "cache"),
    )

    got = lfs.get("foo.grd")

    assert got == str(good_payload)
    assert calls == [
        "https://nadc.example/res/file_upload/download?id=1",
        "https://mirror.example/nlte_grids/foo_v1.grd",
    ]


def test_get_falls_back_when_first_payload_size_mismatches(monkeypatch, tmp_path):
    bad_payload = tmp_path / "bad.grd"
    bad_payload.write_bytes(b"tiny")
    good_payload = tmp_path / "good.grd"
    good_payload.write_bytes(b"grid-data")

    calls = []

    def fake_download_file(url, cache=True, pkgname=""):
        calls.append(url)
        if "nadc" in url:
            return str(bad_payload)
        return str(good_payload)

    monkeypatch.setattr(lfs_module, "download_file", fake_download_file)

    lfs = LargeFileStorage(
        server=["https://mirror.example"],
        pointers={
            "foo.grd": [
                {
                    "url": "https://nadc.example/res/file_upload/download?id=1",
                    "size": 12,
                },
                {"url": "nlte_grids/foo_v1.grd", "size": len(good_payload.read_bytes())},
            ]
        },
        storage=str(tmp_path / "cache"),
    )

    got = lfs.get("foo.grd")

    assert got == str(good_payload)
    assert calls == [
        "https://nadc.example/res/file_upload/download?id=1",
        "https://mirror.example/nlte_grids/foo_v1.grd",
    ]


def test_get_falls_back_when_first_payload_checksum_mismatches(monkeypatch, tmp_path):
    bad_payload = tmp_path / "bad.grd"
    bad_payload.write_bytes(b"grid-a")
    good_payload = tmp_path / "good.grd"
    good_payload.write_bytes(b"grid-b")
    good_sha256 = hashlib.sha256(good_payload.read_bytes()).hexdigest()

    calls = []

    def fake_download_file(url, cache=True, pkgname=""):
        calls.append(url)
        if "nadc" in url:
            return str(bad_payload)
        return str(good_payload)

    monkeypatch.setattr(lfs_module, "download_file", fake_download_file)

    lfs = LargeFileStorage(
        server=["https://mirror.example"],
        pointers={
            "foo.grd": [
                {
                    "url": "https://nadc.example/res/file_upload/download?id=1",
                    "sha256": "deadbeef",
                },
                {
                    "url": "nlte_grids/foo_v1.grd",
                    "sha256": good_sha256,
                },
            ]
        },
        storage=str(tmp_path / "cache"),
    )

    got = lfs.get("foo.grd")

    assert got == str(good_payload)
    assert calls == [
        "https://nadc.example/res/file_upload/download?id=1",
        "https://mirror.example/nlte_grids/foo_v1.grd",
    ]


def test_get_falls_back_when_first_payload_md5_mismatches(monkeypatch, tmp_path):
    bad_payload = tmp_path / "bad.grd"
    bad_payload.write_bytes(b"grid-a")
    good_payload = tmp_path / "good.grd"
    good_payload.write_bytes(b"grid-b")
    good_md5 = hashlib.md5(good_payload.read_bytes()).hexdigest()

    calls = []

    def fake_download_file(url, cache=True, pkgname=""):
        calls.append(url)
        if "zenodo" in url:
            return str(bad_payload)
        return str(good_payload)

    monkeypatch.setattr(lfs_module, "download_file", fake_download_file)

    lfs = LargeFileStorage(
        server=["https://mirror.example"],
        pointers={
            "foo.grd": [
                {
                    "url": "https://zenodo.org/records/1/files/foo.grd?download=1",
                    "md5": "deadbeef",
                },
                {
                    "url": "nlte_grids/foo_v1.grd",
                    "md5": good_md5,
                },
            ]
        },
        storage=str(tmp_path / "cache"),
    )

    got = lfs.get("foo.grd")

    assert got == str(good_payload)
    assert calls == [
        "https://zenodo.org/records/1/files/foo.grd?download=1",
        "https://mirror.example/nlte_grids/foo_v1.grd",
    ]


def test_get_rejects_when_md5_matches_but_sha256_mismatches(monkeypatch, tmp_path):
    payload = tmp_path / "payload.grd"
    payload.write_bytes(b"grid-a")
    good_md5 = hashlib.md5(payload.read_bytes()).hexdigest()

    monkeypatch.setattr(
        lfs_module,
        "download_file",
        lambda url, cache=True, pkgname="": str(payload),
    )

    lfs = LargeFileStorage(
        server=["https://mirror.example"],
        pointers={
            "foo.grd": {
                "url": "https://zenodo.org/records/1/files/foo.grd?download=1",
                "md5": good_md5,
                "sha256": "deadbeef",
            }
        },
        storage=str(tmp_path / "cache"),
    )

    with pytest.raises(FileNotFoundError, match=r"checksum mismatch"):
        lfs.get("foo.grd")


def _make_test_tarball(tmp_path, root_name, members):
    root = tmp_path / root_name
    root.mkdir()
    for name, data in members.items():
        target = root / name
        target.write_bytes(data)

    tar_path = tmp_path / f"{root_name}.tar.gz"
    with tarfile.open(tar_path, "w:gz") as tf:
        tf.add(root, arcname=root_name)
    return tar_path


def test_get_extracts_historical_nlte_tar_archives(monkeypatch, tmp_path):
    tar_path = _make_test_tarball(
        tmp_path,
        "nlte_Na_scatt_pysme",
        {
            "label_Na.txt": b"label",
            "atmos_Na.txt": b"atmos",
            "nlte_Na_scatt_pysme.grd": b"grid-data",
        },
    )

    monkeypatch.setattr(
        lfs_module,
        "download_file",
        lambda url, cache=True, pkgname="": str(tar_path),
    )

    def fake_store_processed_file(self, url, filename):
        cached = tmp_path / "cached_nlte_Na_scatt_pysme.grd"
        shutil.copy2(filename, cached)
        return str(cached)

    monkeypatch.setattr(
        LargeFileStorage,
        "_store_processed_file",
        fake_store_processed_file,
    )

    lfs = LargeFileStorage(
        server=["https://mirror-a.example"],
        pointers={
            "nlte_Na_pysme.grd": "https://zenodo.org/api/records/3960831/files/nlte_Na_scatt_pysme.tar.gz/content"
        },
        storage=str(tmp_path / "cache"),
    )

    got = Path(lfs.get("nlte_Na_pysme.grd"))

    assert got.name == "cached_nlte_Na_scatt_pysme.grd"
    assert got.read_bytes() == b"grid-data"


def test_get_rejects_historical_nlte_tar_without_matching_pysme_grid(monkeypatch, tmp_path):
    tar_path = _make_test_tarball(
        tmp_path,
        "nlte_Na_scatt_pysme",
        {
            "label_Na.txt": b"label",
            "atmos_Na.txt": b"atmos",
        },
    )

    monkeypatch.setattr(
        lfs_module,
        "download_file",
        lambda url, cache=True, pkgname="": str(tar_path),
    )

    lfs = LargeFileStorage(
        server=["https://mirror-a.example"],
        pointers={
            "nlte_Na_pysme.grd": "https://zenodo.org/api/records/3960831/files/nlte_Na_scatt_pysme.tar.gz/content"
        },
        storage=str(tmp_path / "cache"),
    )

    with pytest.raises(FileNotFoundError, match=r"Expected exactly one pysme \.grd"):
        lfs.get("nlte_Na_pysme.grd")


def test_get_rejects_historical_nlte_tar_with_multiple_matching_pysme_grids(monkeypatch, tmp_path):
    tar_path = _make_test_tarball(
        tmp_path,
        "nlte_Na_scatt_pysme",
        {
            "label_Na.txt": b"label",
            "atmos_Na.txt": b"atmos",
            "nlte_Na_scatt_pysme.grd": b"grid-a",
            "nlte_Na_ama51_pysme.grd": b"grid-b",
        },
    )

    monkeypatch.setattr(
        lfs_module,
        "download_file",
        lambda url, cache=True, pkgname="": str(tar_path),
    )

    lfs = LargeFileStorage(
        server=["https://mirror-a.example"],
        pointers={
            "nlte_Na_pysme.grd": "https://zenodo.org/api/records/3960831/files/nlte_Na_scatt_pysme.tar.gz/content"
        },
        storage=str(tmp_path / "cache"),
    )

    with pytest.raises(FileNotFoundError, match=r"Expected exactly one pysme \.grd"):
        lfs.get("nlte_Na_pysme.grd")
