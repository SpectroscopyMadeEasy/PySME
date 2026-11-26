#!/usr/bin/env python3
"""Convert one or more pytest JUnit reports to compact JSON and a notice."""

from __future__ import annotations

import argparse
import json
import re
import xml.etree.ElementTree as ET
from pathlib import Path

from annotations import emit


def testcase_name(case: ET.Element, suite: str) -> str:
    classname = case.attrib.get("classname", "")
    name = case.attrib.get("name", "unknown")
    prefix = f"{suite}::" if suite else ""
    return f"{prefix}{classname}::{name}".strip(":")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--platform", required=True)
    parser.add_argument("reports", type=Path, nargs="+")
    args = parser.parse_args()
    counts = {"passed": 0, "failed": 0, "errors": 0, "skipped": 0}
    failed_tests: list[str] = []
    error_tests: list[str] = []
    missing: list[str] = []
    for report_path in args.reports:
        if not report_path.exists():
            missing.append(str(report_path))
            continue
        root = ET.parse(report_path).getroot()
        suite_label = re.sub(r"-(windows|linux|macos)$", "", report_path.stem)
        for case in root.iter("testcase"):
            if case.find("failure") is not None:
                counts["failed"] += 1
                failed_tests.append(testcase_name(case, suite_label))
            elif case.find("error") is not None:
                counts["errors"] += 1
                error_tests.append(testcase_name(case, suite_label))
            elif case.find("skipped") is not None:
                counts["skipped"] += 1
            else:
                counts["passed"] += 1
    payload = {
        "platform": args.platform,
        "counts": counts,
        "failed_tests": sorted(failed_tests),
        "error_tests": sorted(error_tests),
        "missing_reports": missing,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    emit(
        "notice",
        f"Pytest {args.platform}",
        f"passed={counts['passed']} failed={counts['failed']} errors={counts['errors']} skipped={counts['skipped']} "
        f"failed_tests={(failed_tests + error_tests)[:10] or 'none'} missing={missing or 'none'}",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
