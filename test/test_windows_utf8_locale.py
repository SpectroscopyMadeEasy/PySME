import ctypes
import locale
import logging
from types import SimpleNamespace
from unittest.mock import Mock

import pytest

from pysme import sme_synth


@pytest.fixture(autouse=True)
def reset_warning_state(monkeypatch):
    monkeypatch.setattr(sme_synth, "_WINDOWS_UTF8_LOCALE_WARNING_EMITTED", False)


def mock_ucrt_setlocale(monkeypatch, result):
    setlocale = Mock(return_value=result)
    load_ucrt = Mock(return_value=SimpleNamespace(setlocale=setlocale))
    monkeypatch.setattr(ctypes, "CDLL", load_ucrt)
    return load_ucrt, setlocale


def test_windows_utf8_locale_is_on_by_default(monkeypatch):
    monkeypatch.setattr(sme_synth.sys, "platform", "win32")
    monkeypatch.delenv(sme_synth._WINDOWS_UTF8_LOCALE_ENV, raising=False)
    load_ucrt, setlocale = mock_ucrt_setlocale(monkeypatch, b".UTF-8")

    sme_synth._windows_enable_utf8_locale(r"C:\pysme\data")

    load_ucrt.assert_called_once_with("ucrtbase.dll")
    setlocale.assert_called_once_with(locale.LC_CTYPE, b".UTF-8")


@pytest.mark.parametrize("value", ["0", "false", "NO", " off "])
def test_windows_utf8_locale_can_be_disabled(monkeypatch, value):
    monkeypatch.setattr(sme_synth.sys, "platform", "win32")
    monkeypatch.setenv(sme_synth._WINDOWS_UTF8_LOCALE_ENV, value)
    load_ucrt, _ = mock_ucrt_setlocale(monkeypatch, b".UTF-8")

    sme_synth._windows_enable_utf8_locale(r"C:\pysme\data")

    load_ucrt.assert_not_called()


def test_windows_utf8_locale_failure_warns_once_for_ascii_path(monkeypatch, caplog):
    monkeypatch.setattr(sme_synth.sys, "platform", "win32")
    monkeypatch.delenv(sme_synth._WINDOWS_UTF8_LOCALE_ENV, raising=False)
    _, setlocale = mock_ucrt_setlocale(monkeypatch, None)

    with caplog.at_level(logging.WARNING, logger=sme_synth.__name__):
        sme_synth._windows_enable_utf8_locale(r"C:\pysme\data")
        sme_synth._windows_enable_utf8_locale(r"C:\pysme\other-data")

    assert setlocale.call_count == 2
    warnings = [
        record
        for record in caplog.records
        if record.name == sme_synth.__name__
        and "could not set the Windows UCRT" in record.message
    ]
    assert len(warnings) == 1


def test_windows_utf8_locale_failure_rejects_non_ascii_path(monkeypatch):
    monkeypatch.setattr(sme_synth.sys, "platform", "win32")
    monkeypatch.delenv(sme_synth._WINDOWS_UTF8_LOCALE_ENV, raising=False)
    mock_ucrt_setlocale(monkeypatch, None)

    with pytest.raises(RuntimeError, match="non-ASCII.*Windows 10 version 1803"):
        sme_synth._windows_enable_utf8_locale(r"C:\测试用户\pysme\data")


def test_non_windows_does_not_touch_the_locale(monkeypatch):
    monkeypatch.setattr(sme_synth.sys, "platform", "linux")
    monkeypatch.delenv(sme_synth._WINDOWS_UTF8_LOCALE_ENV, raising=False)
    load_ucrt, _ = mock_ucrt_setlocale(monkeypatch, b".UTF-8")

    sme_synth._windows_enable_utf8_locale("/tmp/café/pysme/data")

    load_ucrt.assert_not_called()
