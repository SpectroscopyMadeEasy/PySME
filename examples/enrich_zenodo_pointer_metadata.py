#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Backfill Zenodo md5/size metadata into PySME pointer files.

This script scans a pointer JSON file, finds Zenodo-backed download targets,
queries the Zenodo Records API for file metadata, and upgrades matching
pointer targets to dict entries with `url`, `md5`, and `size`.

Validation fields refer to the downloadable object itself, e.g. a Zenodo
`.tar.gz` URL is annotated with the tarball's checksum and size.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any
from urllib.error import HTTPError, URLError
from urllib.parse import urlparse, unquote
from urllib.request import urlopen


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "pointer_file",
        nargs="?",
        default=Path(__file__).resolve().parents[1] / "src/pysme/datafiles_nlte.json",
        type=Path,
        help="Pointer JSON to enrich (default: src/pysme/datafiles_nlte.json).",
    )
    parser.add_argument(
        "--write",
        action="store_true",
        help="Write changes back to the pointer file in place.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing md5/size fields instead of leaving them unchanged.",
    )
    return parser.parse_args()


def load_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as f:
        return json.load(f)


def dump_json(path: Path, payload: dict[str, Any]) -> None:
    with path.open("w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2, sort_keys=False)
        f.write("\n")


def is_zenodo_target(target: Any) -> bool:
    if isinstance(target, str):
        return "zenodo.org" in target
    if isinstance(target, dict):
        return "zenodo.org" in str(target.get("url", ""))
    return False


def target_url(target: Any) -> str:
    if isinstance(target, str):
        return target
    if isinstance(target, dict):
        return str(target["url"])
    raise TypeError(f"Unsupported target type: {type(target)!r}")


def parse_zenodo_record_target(url: str) -> tuple[str, str]:
    parsed = urlparse(url)
    parts = [part for part in parsed.path.split("/") if part]

    if len(parts) >= 6 and parts[:3] == ["api", "records", parts[2]] and parts[3] == "files":
        return parts[2], unquote(parts[4])

    if len(parts) >= 4 and parts[0] == "records" and parts[2] == "files":
        return parts[1], unquote(parts[3])

    raise ValueError(f"Unsupported Zenodo URL format: {url}")


def fetch_record_metadata(record_id: str) -> dict[str, Any]:
    url = f"https://zenodo.org/api/records/{record_id}"
    with urlopen(url) as response:
        return json.load(response)


def extract_file_metadata(record_payload: dict[str, Any], filename: str) -> tuple[str | None, int | None]:
    files = record_payload.get("files", [])
    for entry in files:
        key = entry.get("key") or entry.get("filename")
        if key != filename:
            continue

        checksum = entry.get("checksum")
        md5 = None
        if isinstance(checksum, str):
            md5 = checksum.split(":", 1)[1] if checksum.startswith("md5:") else checksum

        size = entry.get("size")
        if size is None and "filesize" in entry:
            try:
                size = int(entry["filesize"])
            except Exception:
                size = entry["filesize"]
        return md5, size

    raise KeyError(f"Could not find file {filename!r} in Zenodo record response")


def upgrade_target(target: Any, md5: str | None, size: int | None, overwrite: bool) -> dict[str, Any]:
    base = {"url": target_url(target)}
    if isinstance(target, dict):
        base.update(target)

    if md5 is not None and (overwrite or "md5" not in base):
        base["md5"] = md5
    if size is not None and (overwrite or "size" not in base):
        base["size"] = size
    return base


def enrich_pointer_data(data: dict[str, Any], overwrite: bool) -> tuple[dict[str, Any], list[str]]:
    updated = {}
    log_lines: list[str] = []
    record_cache: dict[str, dict[str, Any]] = {}

    for key, value in data.items():
        items = value if isinstance(value, list) else [value]
        new_items = []

        for item in items:
            if not is_zenodo_target(item):
                new_items.append(item)
                continue

            url = target_url(item)
            record_id, filename = parse_zenodo_record_target(url)
            if record_id not in record_cache:
                record_cache[record_id] = fetch_record_metadata(record_id)
            md5, size = extract_file_metadata(record_cache[record_id], filename)
            enriched = upgrade_target(item, md5, size, overwrite=overwrite)
            new_items.append(enriched)
            log_lines.append(f"{key}: {filename} (record {record_id}) md5={md5} size={size}")

        updated[key] = new_items if isinstance(value, list) else new_items[0]

    return updated, log_lines


def main() -> int:
    args = parse_args()
    pointer_file = args.pointer_file.resolve(strict=False)
    data = load_json(pointer_file)

    try:
        updated, log_lines = enrich_pointer_data(data, overwrite=args.overwrite)
    except (HTTPError, URLError, ValueError, KeyError) as exc:
        print(f"[error] {exc}", file=sys.stderr)
        return 1

    if not log_lines:
        print("[info] no Zenodo targets found")
        return 0

    for line in log_lines:
        print(f"[match] {line}")

    if args.write:
        dump_json(pointer_file, updated)
        print(f"[done] wrote {pointer_file}")
    else:
        print("[done] dry-run only; re-run with --write to update the pointer file")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
