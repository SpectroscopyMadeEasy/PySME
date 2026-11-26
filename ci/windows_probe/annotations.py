#!/usr/bin/env python3
"""Emit compact, correctly escaped GitHub workflow annotations."""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path


MAX_MESSAGE_CHARS = 3800
MAX_TAIL_CHUNKS = 3
ERROR_CONTEXT_BEFORE = 2
ERROR_CONTEXT_AFTER = 4

# Match actual diagnostic forms, not arbitrary occurrences of "error" in flag
# names such as -Werror or -Wformat-truncation.
REAL_ERROR_PATTERNS = (
    re.compile(r"^\s*(?:cmake\s+)?error(?:\s+at\b|\s+[A-Z]+\d+\s*:|\s*:)", re.I),
    re.compile(r"^.+:\d+(?::\d+)?:\s*(?:fatal\s+)?error(?:\s+[A-Z]+\d+)?\s*:", re.I),
    re.compile(r"^.+\(\d+(?:,\d+)?\)\s*:\s*(?:fatal\s+)?error\s+[A-Z]+\d+\s*:", re.I),
    re.compile(r"^\s*[A-Z0-9_.+-]+\s*:\s*(?:fatal\s+)?error(?:\s+[A-Z]+\d+)?\s*:", re.I),
)


def escape(value: object) -> str:
    return str(value).replace("%", "%25").replace("\r", "%0D").replace("\n", "%0A")


def emit(level: str, title: str, message: object) -> None:
    line = f"::{level} title={escape(title)}::{escape(message)}\n"
    stream = getattr(sys.stdout, "buffer", None)
    if stream is None:
        # StringIO and similar test streams have no byte interface. A real
        # console takes the byte path, bypassing its possibly narrow codec.
        sys.stdout.write(line)
        sys.stdout.flush()
    else:
        stream.write(line.encode("utf-8"))
        stream.flush()


def is_real_error(line: str) -> bool:
    return any(pattern.search(line) for pattern in REAL_ERROR_PATTERNS)


def split_from_end(message: str, limit: int, max_chunks: int) -> list[str]:
    """Return newest text first, with each escaped chunk within limit."""
    chunks: list[str] = []
    end = len(message)
    while end and len(chunks) < max_chunks:
        low, high = 0, end
        while low < high:
            middle = (low + high) // 2
            if len(escape(message[middle:end])) <= limit:
                high = middle
            else:
                low = middle + 1
        start = low
        chunks.append(message[start:end].strip("\n") or "(blank log tail)")
        end = start
    return chunks


def truncate_from_start(message: str, limit: int) -> str:
    if len(escape(message)) <= limit:
        return message
    low, high = 0, len(message)
    while low < high:
        middle = (low + high + 1) // 2
        if len(escape(message[:middle])) <= limit:
            low = middle
        else:
            high = middle - 1
    return message[:low]


def tail_chunks(lines: list[str], tail: int) -> list[str]:
    if not lines:
        return ["(log was empty)"]
    return split_from_end("\n".join(lines[-tail:]), MAX_MESSAGE_CHARS, MAX_TAIL_CHUNKS)


def first_error_context(lines: list[str]) -> tuple[int, str] | None:
    for index, line in enumerate(lines):
        if is_real_error(line):
            start = max(0, index - ERROR_CONTEXT_BEFORE)
            end = min(len(lines), index + ERROR_CONTEXT_AFTER + 1)
            before = lines[start:index]
            after = lines[index + 1 : end]
            context = "\n".join(
                [
                    f"first real error: {line}",
                    *("context before: " + value for value in before),
                    *("context after: " + value for value in after),
                ]
            )
            return index, truncate_from_start(context, MAX_MESSAGE_CHARS)
    return None


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
    lines = args.log.read_text(encoding="utf-8", errors="replace").splitlines()
    chunks = tail_chunks(lines, args.tail)
    for index, chunk in enumerate(chunks, start=1):
        suffix = " (newest)" if index == 1 else ""
        emit(args.level, f"{args.title} tail {index}/{len(chunks)}{suffix}", chunk)

    error = first_error_context(lines)
    tail_start = max(0, len(lines) - args.tail)
    if error is not None and error[0] < tail_start:
        emit(args.level, f"{args.title} first error", error[1])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
