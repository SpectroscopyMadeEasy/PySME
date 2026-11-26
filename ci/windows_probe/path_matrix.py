#!/usr/bin/env python3
"""Run independent install/profile path combinations against an ASCII result."""

from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

from annotations import emit


CASES = ("solar_fixed", "hydrogen_fixed")
SHORT_PATH_ENV = "PYSME_SMELIB_WINDOWS_SHORT_PATH"
UTF8_LOCALE_ENV = "PYSME_SMELIB_WINDOWS_UTF8_LOCALE"
FAILED_OPEN_RE = re.compile(
    r"(?:tables?\s+not\s+found|did\s+not\s+find[^\r\n]*\bfile\b|"
    r"(?:cannot|can't|could\s+not|unable\s+to|failed\s+to)\s+open|"
    r"error\s+(?:reading|opening|open)[^\r\n]*\bfile\b|"
    r"eof\s+found\s+while\s+reading[^\r\n]*\bfile\b)",
    re.IGNORECASE,
)
PROFILE_PATHS = {
    "data.hlineprof": "~/.sme/hlineprof",
    "data.atmospheres": "~/.sme/atmospheres",
    "data.nlte_grids": "~/.sme/nlte_grids",
    "data.cache.atmospheres": "~/.sme/atmospheres/cache",
    "data.cache.nlte_grids": "~/.sme/nlte_grids/cache",
}


def compact(value: object, limit: int = 260) -> str:
    text = " ".join(str(value).split())
    return text if len(text) <= limit else text[: limit - 3] + "..."


def current_ucrt_lc_ctype() -> str | None:
    if sys.platform != "win32":
        return None
    import ctypes
    import locale

    ucrt = ctypes.CDLL("ucrtbase.dll")
    setlocale = ucrt.setlocale
    setlocale.argtypes = (ctypes.c_int, ctypes.c_char_p)
    setlocale.restype = ctypes.c_char_p
    value = setlocale(locale.LC_CTYPE, None)
    return None if value is None else value.decode("ascii", errors="replace")


def copy_profile(source_profile: Path, target_profile: Path) -> None:
    source = source_profile.resolve()
    target = target_profile.resolve()
    source_sme = source / ".sme"
    target_sme = target / ".sme"
    if not source_sme.is_dir():
        raise FileNotFoundError(f"prefetched profile is missing {source_sme}")
    target.mkdir(parents=True, exist_ok=True)
    if target_sme.exists():
        shutil.rmtree(target_sme)
    shutil.copytree(source_sme, target_sme)

    config_path = target_sme / "config.json"
    config = json.loads(config_path.read_text(encoding="utf-8"))
    config.update(PROFILE_PATHS)
    config_path.write_text(
        json.dumps(config, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )

    # A copied cache must not retain metadata that points back to the source
    # profile. Check UTF-8 and UTF-16LE because either may occur on Windows.
    needles = (str(source).encode("utf-8"), str(source).encode("utf-16le"))
    offenders = []
    for entry in target_sme.rglob("*"):
        if entry.is_file():
            data = entry.read_bytes()
            if any(needle and needle in data for needle in needles):
                offenders.append(str(entry.relative_to(target_sme)))
    if offenders:
        raise RuntimeError(
            "copied cache contains absolute source-profile paths: "
            + ", ".join(offenders[:5])
        )


def case_main(args: argparse.Namespace) -> int:
    result = {"case": args.case}
    try:
        import numpy as np
        from pysme.config import Config
        from pysme.smelib.libtools import get_full_datadir
        from reference_calculation import synthesize_case

        profile = Path.home().resolve()
        config = Config()
        for key in PROFILE_PATHS:
            resolved = Path(config[key]).expanduser().resolve()
            if not resolved.is_relative_to(profile):
                raise RuntimeError(
                    f"{key} escaped row profile: {resolved} is not below {profile}"
                )
        datadir = Path(get_full_datadir()).resolve()
        if not datadir.is_relative_to(Path(sys.prefix).resolve()):
            raise RuntimeError(
                f"SMElib data directory escaped row install: {datadir} is not below {sys.prefix}"
            )
        with np.load(args.baseline, allow_pickle=False) as baseline:
            baseline_wave = np.asarray(baseline[f"{args.case}__wave"], dtype=float)
            baseline_flux = np.asarray(baseline[f"{args.case}__flux"], dtype=float)
        wave, flux, _ = synthesize_case(args.repo.resolve(), args.case)
        if not np.array_equal(wave, baseline_wave):
            raise ValueError("wavelength grid differs from the ASCII baseline")
        difference = np.asarray(flux - baseline_flux, dtype=float)
        if not np.all(np.isfinite(difference)):
            raise RuntimeError("flux difference contains non-finite values")
        max_abs = float(np.max(np.abs(difference)))
        rms = float(np.sqrt(np.mean(difference**2)))
        if max_abs > 1e-12:
            raise ValueError(f"max |dflux|={max_abs:.9g} exceeds 1e-12")
        result.update(
            status="PASS",
            max_abs=max_abs,
            rms=rms,
            ucrt_lc_ctype=current_ucrt_lc_ctype(),
        )
    except BaseException as exc:
        result.update(status="FAIL", exception=type(exc).__name__, message=str(exc))
    args.result.parent.mkdir(parents=True, exist_ok=True)
    args.result.write_text(json.dumps(result, ensure_ascii=False), encoding="utf-8")
    return 0 if result["status"] == "PASS" else 1


def output_tail(completed: subprocess.CompletedProcess[str], lines: int = 4) -> str:
    combined = "\n".join(part for part in (completed.stdout, completed.stderr) if part)
    values = [line.strip() for line in combined.splitlines() if line.strip()]
    return " | ".join(values[-lines:]) or "(no process output)"


def row_main(args: argparse.Namespace) -> int:
    args.work.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    env["HOME"] = str(args.profile)
    env["USERPROFILE"] = str(args.profile)
    env["PYTHONIOENCODING"] = "utf-8"
    env["GFORTRAN_UNBUFFERED_ALL"] = "y"
    if args.short_path:
        env[SHORT_PATH_ENV] = "1"
    else:
        env.pop(SHORT_PATH_ENV, None)
    if args.utf8_locale:
        env[UTF8_LOCALE_ENV] = "1"
    else:
        env.pop(UTF8_LOCALE_ENV, None)

    outcomes = []
    for case in CASES:
        result_path = args.work / f"{args.label}-{case}.json"
        log_path = args.work / f"{args.label}-{case}.log"
        command = [
            sys.executable,
            "-u",
            str(Path(__file__).resolve()),
            "--case",
            case,
            "--baseline",
            str(args.baseline),
            "--repo",
            str(args.repo),
            "--result",
            str(result_path),
        ]
        completed = subprocess.run(
            command, env=env, capture_output=True, text=True,
            encoding="utf-8", errors="replace"
        )
        log_path.write_text(
            completed.stdout + completed.stderr, encoding="utf-8", errors="replace"
        )
        captured = completed.stdout + completed.stderr
        if result_path.exists():
            result = json.loads(result_path.read_text(encoding="utf-8"))
        else:
            result = {
                "case": case,
                "status": "FAIL",
                "exception": "ProcessExit",
                "message": f"return code {completed.returncode}",
            }
        failed_open = FAILED_OPEN_RE.search(captured)
        if failed_open is not None:
            result.update(
                status="FAIL",
                exception="NativeFileOpenFailure",
                message=f"captured native failed-open diagnostic: {failed_open.group(0)}",
            )
        result["tail"] = output_tail(completed)
        outcomes.append(result)

    failed = [value for value in outcomes if value["status"] != "PASS"]
    if not failed:
        detail = "; ".join(
            f"{value['case']} max|dflux|={value['max_abs']:.9g}" for value in outcomes
        )
        locales = {value.get("ucrt_lc_ctype") for value in outcomes}
        if locales != {None}:
            detail += f"; UCRT_LC_CTYPE={compact(sorted(locales, key=str))}"
        emit("notice", f"Path matrix {args.label} PASS", detail)
        return 0

    parts = []
    for value in outcomes:
        if value["status"] == "PASS":
            parts.append(f"{value['case']} PASS max={value['max_abs']:.3g}")
        else:
            parts.append(
                f"{value['case']} FAIL {value['exception']}: {compact(value['message'])}; "
                f"tail={compact(value['tail'])}"
            )
    emit("notice", f"Path matrix {args.label} FAIL", "; ".join(parts))
    return 1


def controller_main(args: argparse.Namespace) -> int:
    copy_profile(args.source_profile, args.ascii_profile)
    copy_profile(args.source_profile, args.chinese_profile)
    rows = (
        (
            "a-install-ascii-profile-ascii",
            args.ascii_python, args.ascii_profile, False, False,
        ),
        (
            "b-install-ascii-profile-chinese",
            args.ascii_python, args.chinese_profile, False, False,
        ),
        (
            "c-install-chinese-profile-ascii",
            args.chinese_python, args.ascii_profile, False, False,
        ),
        (
            "d-both-chinese",
            args.chinese_python, args.chinese_profile, False, False,
        ),
        (
            "d-both-chinese-short-path",
            args.chinese_python, args.chinese_profile, True, False,
        ),
        (
            "d-both-chinese-utf8-locale",
            args.chinese_python, args.chinese_profile, False, True,
        ),
    )
    failed = False
    for label, python, profile, short_path, utf8_locale in rows:
        command = [
            str(python), "-u", str(Path(__file__).resolve()), "--row",
            "--label", label, "--profile", str(profile),
            "--baseline", str(args.baseline), "--repo", str(args.repo),
            "--work", str(args.work),
        ]
        if short_path:
            command.append("--short-path")
        if utf8_locale:
            command.append("--utf8-locale")
        completed = subprocess.run(
            command, capture_output=True, text=True, encoding="utf-8", errors="replace"
        )
        if completed.stdout:
            sys.stdout.buffer.write(completed.stdout.encode("utf-8"))
            sys.stdout.buffer.flush()
        if completed.returncode != 0:
            failed = True
            if "::notice " not in completed.stdout:
                emit(
                    "notice", f"Path matrix {label} FAIL",
                    f"ProcessExit: row process returned {completed.returncode}; "
                    f"tail={compact(output_tail(completed))}",
                )
    return int(failed)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--row", action="store_true")
    mode.add_argument("--case", choices=CASES)
    parser.add_argument("--label")
    parser.add_argument("--profile", type=Path)
    parser.add_argument("--baseline", type=Path, required=True)
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--work", type=Path)
    parser.add_argument("--result", type=Path)
    parser.add_argument("--short-path", action="store_true")
    parser.add_argument("--utf8-locale", action="store_true")
    parser.add_argument("--source-profile", type=Path)
    parser.add_argument("--ascii-profile", type=Path)
    parser.add_argument("--chinese-profile", type=Path)
    parser.add_argument("--ascii-python", type=Path)
    parser.add_argument("--chinese-python", type=Path)
    args = parser.parse_args()
    if args.case is not None:
        required = (args.result,)
    elif args.row:
        required = (args.label, args.profile, args.work)
    else:
        required = (
            args.source_profile, args.ascii_profile, args.chinese_profile,
            args.ascii_python, args.chinese_python, args.work,
        )
    if any(value is None for value in required):
        parser.error("missing arguments required by the selected mode")
    return args


def main() -> int:
    args = parse_args()
    if args.case is not None:
        return case_main(args)
    if args.row:
        return row_main(args)
    return controller_main(args)


if __name__ == "__main__":
    raise SystemExit(main())
