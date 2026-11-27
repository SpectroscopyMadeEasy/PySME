# -*- coding: utf-8 -*-
"""The fallback SMElib download must be the release of the bundled smelib submodule."""
import subprocess
from pathlib import Path

import pytest

from pysme.smelib import libtools

ROOT = Path(__file__).resolve().parents[1]


def test_fallback_release_matches_the_smelib_submodule():
    try:
        tag = subprocess.run(
            ["git", "-C", str(ROOT / "smelib"), "describe", "--tags", "--exact-match"],
            capture_output=True, text=True, check=True,
        ).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        pytest.skip("smelib submodule is not a git checkout at a release tag")
    assert libtools.SMELIB_RELEASE == tag
    assert libtools.smelib_releases["default"] == f"download/{tag}"
