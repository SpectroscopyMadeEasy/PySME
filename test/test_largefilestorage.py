# -*- coding: utf-8 -*-
import pytest
import requests
from astropy.utils.data import is_url_in_cache

from pysme.large_file_storage import setup_atmo, setup_nlte


def lfs_available():
    lfs_atmo = setup_atmo()
    urls = lfs_atmo.get_urls("marcs2012.sav")

    if any(is_url_in_cache(url, pkgname=lfs_atmo.PKGNAME) for url in urls):
        return True

    for url in urls:
        try:
            with requests.get(
                url,
                headers={"Range": "bytes=0-1"},
                stream=True,
                timeout=5,
            ) as response:
                if response.status_code in (200, 206):
                    if response.raw.read(2) == b"\x1f\x8b":
                        return True
        except requests.RequestException:
            continue

    return False


skipif_lfs = pytest.mark.skipif(not lfs_available(), reason="LFS not available")


@pytest.fixture
def lfs_nlte():
    lfs_nlte = setup_nlte()
    yield lfs_nlte


@pytest.fixture
def lfs_atmo():
    lfs_atmo = setup_atmo()
    yield lfs_atmo
