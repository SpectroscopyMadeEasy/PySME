#!/usr/bin/env python3
"""Download all tracked PySME large data files (atmospheres/NLTE)."""

import argparse
import json
from pathlib import Path

from pysme.config import Config
from pysme.large_file_storage import setup_atmo, setup_nlte


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--kind",
        choices=["all", "atmo", "nlte"],
        default="all",
        help="Which data groups to download",
    )
    parser.add_argument(
        "--config",
        default="~/.sme/config.json",
        help="Path to PySME config JSON",
    )
    parser.add_argument(
        "--manifest",
        default="download_large_datafiles_manifest.json",
        help="Output manifest path",
    )
    parser.add_argument(
        "--fail-fast",
        action="store_true",
        help="Stop on first failed file instead of continuing",
    )
    return parser.parse_args()


def fetch_group(group_name, lfs, fail_fast):
    results = []
    for key in sorted(lfs.pointers.keys()):
        item = {"group": group_name, "key": key}
        try:
            path = lfs.get(key)
            item["status"] = "ok"
            item["path"] = path
        except Exception as exc:
            item["status"] = "error"
            item["error"] = str(exc)
            if fail_fast:
                results.append(item)
                return results, False
        results.append(item)
    return results, True


def main():
    args = parse_args()

    cfg = Config(args.config)
    cfg.load()

    groups = []
    if args.kind in ("all", "atmo"):
        groups.append(("atmo", setup_atmo(cfg)))
    if args.kind in ("all", "nlte"):
        groups.append(("nlte", setup_nlte(cfg)))

    all_results = []
    ok = True

    for group_name, lfs in groups:
        print(f"[group] {group_name}: {len(lfs.pointers)} files")
        results, group_ok = fetch_group(group_name, lfs, args.fail_fast)
        all_results.extend(results)
        ok = ok and group_ok
        n_ok = sum(1 for r in results if r["status"] == "ok")
        n_err = sum(1 for r in results if r["status"] == "error")
        print(f"[group] {group_name}: ok={n_ok}, error={n_err}")
        if args.fail_fast and not group_ok:
            break

    manifest_path = Path(args.manifest)
    manifest_path.write_text(json.dumps(all_results, indent=2, ensure_ascii=False) + "\n")
    print(f"[done] Wrote manifest: {manifest_path}")

    if not ok:
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
