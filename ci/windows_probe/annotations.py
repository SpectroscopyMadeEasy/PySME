#!/usr/bin/env python3
"""Emit compact, correctly escaped GitHub workflow annotations."""

from __future__ import annotations

import argparse
import re
from pathlib import Path


def escape(value: object) -> str:
    return str(value).replace("%", "%25").replace("\r", "%0D").replace("\n", "%0A")


def emit(level: str, title: str, message: object) -> None:
    print(f"::{level} title={escape(title)}::{escape(message)}")


def summarize_log(path: Path, tail: int) -> str:
    text = path.read_text(encoding="utf-8", errors="replace")
    lines = text.splitlines()
    first_error = next(
        (
            line
            for line in lines
            if re.search(r"(^|[ :] )(fatal error|error:|error [A-Z0-9]+:)", line, re.I)
        ),
        None,
    )
    selected = lines[-tail:]
    if first_error and first_error not in selected:
        selected = [f"first compiler/error line: {first_error}", *selected]
    return "\n".join(selected) if selected else "(log was empty)"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("log", type=Path)
    parser.add_argument("--title", required=True)
    parser.add_argument("--level", choices=("notice", "warning", "error"), default="error")
    parser.add_argument("--tail", type=int, default=80)
    args = parser.parse_args()
    if not args.log.exists():
        emit(args.level, args.title, f"log not found: {args.log}")
        return 0
    emit(args.level, args.title, summarize_log(args.log, args.tail))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
