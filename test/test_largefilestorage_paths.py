# -*- coding: utf-8 -*-
import os
from pathlib import Path

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
