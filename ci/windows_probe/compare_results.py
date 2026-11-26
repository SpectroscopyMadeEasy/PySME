#!/usr/bin/env python3
"""Compare Windows and macOS probe artifacts against Linux."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from annotations import emit


def load(path: Path) -> dict[str, object]:
    with np.load(path, allow_pickle=False) as data:
        return {key: data[key] for key in data.files}


def scalar_text(value: object) -> str:
    return str(np.asarray(value).item())


def compare_case(candidate: dict[str, object], linux: dict[str, object], name: str):
    wave_ref = np.asarray(linux[f"{name}__wave"], dtype=float)
    flux_ref = np.asarray(linux[f"{name}__flux"], dtype=float)
    wave = np.asarray(candidate[f"{name}__wave"], dtype=float)
    flux = np.asarray(candidate[f"{name}__flux"], dtype=float)
    kind = scalar_text(linux[f"{name}__grid_kind"])
    if kind == "fixed":
        if wave.shape != wave_ref.shape or not np.array_equal(wave, wave_ref):
            raise ValueError(f"{name}: fixed wavelength grids differ")
        difference = flux - flux_ref
        compared = wave.size
        mode = "direct"
    else:
        lo = max(wave[0], wave_ref[0])
        hi = min(wave[-1], wave_ref[-1])
        use = (wave_ref >= lo) & (wave_ref <= hi)
        difference = np.interp(wave_ref[use], wave, flux) - flux_ref[use]
        compared = int(use.sum())
        mode = "interpolated-to-linux"
    return {
        "max": float(np.max(np.abs(difference))),
        "rms": float(np.sqrt(np.mean(difference**2))),
        "candidate_pixels": int(wave.size),
        "linux_pixels": int(wave_ref.size),
        "compared_pixels": compared,
        "kind": kind,
        "mode": mode,
    }


def test_failures(path: Path | None) -> tuple[dict[str, object], set[str]]:
    if path is None or not path.exists():
        return {}, set()
    report = json.loads(path.read_text(encoding="utf-8"))
    failures = set(report.get("failed_tests", [])) | set(report.get("error_tests", []))
    return report, failures


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--windows", type=Path, required=True)
    parser.add_argument("--linux", type=Path, required=True)
    parser.add_argument("--macos", type=Path, required=True)
    parser.add_argument("--windows-tests", type=Path)
    parser.add_argument("--linux-tests", type=Path)
    args = parser.parse_args()
    missing = [str(p) for p in (args.windows, args.linux, args.macos) if not p.exists()]
    if missing:
        emit("error", "Missing probe artifacts", ", ".join(missing))
        return 1

    results = {"windows": load(args.windows), "linux": load(args.linux), "macos": load(args.macos)}
    names = [str(value) for value in np.asarray(results["linux"]["case_names"]).tolist()]
    failures: list[str] = []
    metrics: dict[str, dict[str, dict[str, object]]] = {}
    for name in names:
        metrics[name] = {}
        for platform_name in ("windows", "macos"):
            try:
                metric = compare_case(results[platform_name], results["linux"], name)
            except Exception as exc:
                emit("error", f"Compare {name} {platform_name}-linux", repr(exc))
                failures.append(f"{name}/{platform_name}: {exc}")
                continue
            metrics[name][platform_name] = metric
            emit(
                "notice",
                f"{name} {platform_name}-linux",
                f"grid={metric['kind']} mode={metric['mode']} max={metric['max']:.9g} rms={metric['rms']:.9g} "
                f"pixels={metric['candidate_pixels']}/{metric['linux_pixels']} compared={metric['compared_pixels']}",
            )
        if "windows" in metrics[name] and "macos" in metrics[name] and metrics[name]["windows"]["kind"] == "fixed":
            win_max = float(metrics[name]["windows"]["max"])
            mac_max = float(metrics[name]["macos"]["max"])
            if win_max > 1e-4 and win_max > 2.0 * mac_max:
                message = f"Windows-Linux max {win_max:.9g} > 1e-4 and > 2 x macOS-Linux {mac_max:.9g}"
                emit("error", f"Fixed-grid threshold {name}", message)
                failures.append(f"{name}: {message}")

    fit_messages = []
    for platform_name in ("windows", "linux", "macos"):
        fit = json.loads(scalar_text(results[platform_name]["fit_json"]))
        if fit is None:
            emit("error", f"Fit {platform_name}", "fit result missing")
            failures.append(f"fit/{platform_name}: missing")
        else:
            fit_messages.append(
                f"{platform_name}: value={fit['value']:.9g} uncertainty={fit['uncertainty']:.9g} "
                f"chisq={fit['chisq']:.9g} iterations={fit['iterations']}"
            )
    if fit_messages:
        emit("notice", "Fit comparison", "; ".join(fit_messages))

    win_report, win_failed = test_failures(args.windows_tests)
    linux_report, linux_failed = test_failures(args.linux_tests)
    for platform_name, report in (("windows", win_report), ("linux", linux_report)):
        if not report or report.get("missing_reports"):
            message = f"report={report or 'missing'}"
            emit("error", f"Pytest report incomplete {platform_name}", message)
            failures.append(f"pytest/{platform_name}: {message}")
    windows_only = sorted(win_failed - linux_failed)
    emit(
        "notice",
        "Pytest cross-platform",
        f"windows={win_report.get('counts', 'missing')} linux={linux_report.get('counts', 'missing')} windows_only={windows_only or 'none'}",
    )
    if windows_only:
        failures.extend(f"Windows-only pytest: {name}" for name in windows_only)
        emit("error", "Windows-specific pytest failures", ", ".join(windows_only))
    if failures:
        print("\n".join(failures))
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
