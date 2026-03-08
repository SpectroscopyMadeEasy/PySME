# -*- coding: utf-8 -*-
import os
from pathlib import Path

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
