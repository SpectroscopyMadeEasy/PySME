#!/usr/bin/env python3
"""Load the installed native library and report its resolved Windows modules."""

from __future__ import annotations

import ctypes
import os
import platform
import sys
from pathlib import Path

from annotations import emit


def windows_modules() -> list[str]:
    if os.name != "nt":
        return []
    psapi = ctypes.WinDLL("psapi", use_last_error=True)
    kernel32 = ctypes.WinDLL("kernel32", use_last_error=True)
    kernel32.GetCurrentProcess.restype = ctypes.c_void_p
    psapi.EnumProcessModules.argtypes = (
        ctypes.c_void_p,
        ctypes.POINTER(ctypes.c_void_p),
        ctypes.c_ulong,
        ctypes.POINTER(ctypes.c_ulong),
    )
    psapi.EnumProcessModules.restype = ctypes.c_int
    psapi.GetModuleFileNameExW.argtypes = (
        ctypes.c_void_p,
        ctypes.c_void_p,
        ctypes.POINTER(ctypes.c_wchar),
        ctypes.c_ulong,
    )
    psapi.GetModuleFileNameExW.restype = ctypes.c_ulong
    process = kernel32.GetCurrentProcess()
    modules = (ctypes.c_void_p * 2048)()
    needed = ctypes.c_ulong()
    if not psapi.EnumProcessModules(
        process, modules, ctypes.sizeof(modules), ctypes.byref(needed)
    ):
        raise ctypes.WinError(ctypes.get_last_error())
    count = min(needed.value // ctypes.sizeof(ctypes.c_void_p), len(modules))
    result = []
    for module in modules[:count]:
        buffer = ctypes.create_unicode_buffer(32768)
        if psapi.GetModuleFileNameExW(process, module, buffer, len(buffer)):
            name = buffer.value
            lowered = Path(name).name.lower()
            if any(token in lowered for token in ("sme", "gfortran", "gcc_s", "quadmath", "winpthread")):
                result.append(name)
    return sorted(set(result))


def main() -> int:
    from pysme.sme_synth import ensure_smelib_ready
    from pysme.smelib.libtools import get_full_datadir, get_full_libfile

    native = ensure_smelib_ready()
    print(f"_smelib.LibraryVersion()={native.LibraryVersion()}")
    print(f"SME library={get_full_libfile()}")
    print(f"SME data={get_full_datadir()}")
    modules = windows_modules()
    for module in modules:
        print(f"loaded module={module}")
    emit(
        "notice",
        "Native load",
        f"python={sys.version.split()[0]} arch={platform.machine()} version={native.LibraryVersion()} "
        f"library={get_full_libfile()} loaded_modules={modules or ['reported only on Windows']}",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
