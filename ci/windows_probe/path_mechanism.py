#!/usr/bin/env python3
"""Probe how the Windows UCRT opens ASCII, ACP, and Unicode paths."""

from __future__ import annotations

import ctypes
import locale
import os
import platform
import sys
import tempfile
from pathlib import Path

from annotations import emit


PAYLOAD = b"pysme-ucrt-path-probe\n"


def compact(value: object) -> str:
    return " ".join(str(value).split())


def observed(function) -> object:
    try:
        return function()
    except Exception as exc:
        return f"ERROR({type(exc).__name__}: {exc})"


def python_state(sample: Path) -> dict[str, object]:
    value = 1234.5
    return {
        "preferred": observed(lambda: locale.getpreferredencoding(False)),
        "locale": observed(locale.getlocale),
        "float": observed(lambda: float("1234.5")),
        "format": observed(lambda: format(value, ".2f")),
        "path_roundtrip": observed(lambda: str(Path(str(sample))) == str(sample)),
    }


def open_and_read(open_file, path, mode) -> str:
    stream = open_file(path, mode)
    if not stream:
        return "FAIL(NULL)"
    try:
        buffer = ctypes.create_string_buffer(len(PAYLOAD))
        count = UCRT.fread(buffer, 1, len(PAYLOAD), stream)
        if count != len(PAYLOAD):
            return f"FAIL(read={count}/{len(PAYLOAD)})"
        if buffer.raw != PAYLOAD:
            return "FAIL(content)"
        return "PASS"
    finally:
        UCRT.fclose(stream)


def acp_bytes(path: Path, code_page: int) -> bytes | None:
    try:
        return str(path).encode(f"cp{code_page}", errors="strict")
    except UnicodeEncodeError:
        return None


if sys.platform != "win32":
    emit("notice", "UCRT path mechanism skipped", f"platform={sys.platform}")
    raise SystemExit(0)


KERNEL32 = ctypes.WinDLL("kernel32", use_last_error=True)
KERNEL32.GetACP.argtypes = ()
KERNEL32.GetACP.restype = ctypes.c_uint
KERNEL32.GetOEMCP.argtypes = ()
KERNEL32.GetOEMCP.restype = ctypes.c_uint

UCRT = ctypes.CDLL("ucrtbase.dll")
UCRT.fopen.argtypes = (ctypes.c_char_p, ctypes.c_char_p)
UCRT.fopen.restype = ctypes.c_void_p
UCRT._wfopen.argtypes = (ctypes.c_wchar_p, ctypes.c_wchar_p)
UCRT._wfopen.restype = ctypes.c_void_p
UCRT.fread.argtypes = (
    ctypes.c_void_p,
    ctypes.c_size_t,
    ctypes.c_size_t,
    ctypes.c_void_p,
)
UCRT.fread.restype = ctypes.c_size_t
UCRT.fclose.argtypes = (ctypes.c_void_p,)
UCRT.fclose.restype = ctypes.c_int
UCRT.setlocale.argtypes = (ctypes.c_int, ctypes.c_char_p)
UCRT.setlocale.restype = ctypes.c_char_p

acp = int(KERNEL32.GetACP())
oemcp = int(KERNEL32.GetOEMCP())
root = Path(tempfile.mkdtemp(prefix="pysme-ucrt-paths-"))
directories = (("ascii", "ascii"), ("cp1252", "café"), ("chinese", "测试目录"))
paths: list[tuple[str, Path]] = []
for label, dirname in directories:
    directory = root / dirname
    directory.mkdir()
    path = directory / "probe.bin"
    path.write_bytes(PAYLOAD)
    paths.append((label, path))

baseline_locale = UCRT.setlocale(locale.LC_CTYPE, None)
before = python_state(paths[-1][1])
emit(
    "notice",
    "UCRT path mechanism system",
    f"Windows={platform.platform()} version={platform.version()} "
    f"getwindowsversion={sys.getwindowsversion()} "
    f"GetACP={acp} GetOEMCP={oemcp} initial_LC_CTYPE={baseline_locale!r}",
)

before_results: dict[str, tuple[str, str, str]] = {}
for label, path in paths:
    utf8 = open_and_read(UCRT.fopen, str(path).encode("utf-8"), b"rb")
    encoded = acp_bytes(path, acp)
    ansi = (
        "not representable"
        if encoded is None
        else open_and_read(UCRT.fopen, encoded, b"rb")
    )
    wide = open_and_read(UCRT._wfopen, str(path), "rb")
    before_results[label] = (utf8, ansi, wide)

locale_result = UCRT.setlocale(locale.LC_CTYPE, b".UTF-8")
after = python_state(paths[-1][1])
effective_locale = UCRT.setlocale(locale.LC_CTYPE, None)
changed = {key: before[key] != after[key] for key in before}
emit(
    "notice",
    "UCRT UTF-8 locale",
    f"setlocale={locale_result!r} effective={effective_locale!r} "
    f"changed={compact(changed)} before={compact(before)} after={compact(after)}",
)

control_failed = False
for label, path in paths:
    utf8, ansi, wide = before_results[label]
    post_locale = open_and_read(UCRT.fopen, str(path).encode("utf-8"), b"rb")
    emit(
        "notice",
        f"UCRT path {label}",
        f"utf8-narrow={utf8}; acp-narrow={ansi}; wide={wide}; "
        f"utf8-after-locale={post_locale}",
    )
    control_failed |= wide != "PASS"

if control_failed:
    raise RuntimeError("one or more _wfopen control reads failed")
