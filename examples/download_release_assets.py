#!/usr/bin/env python3
"""Download all assets from a GitHub release tag."""

import argparse
import json
import os
import sys
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen


API_TEMPLATE = "https://api.github.com/repos/{repo}/releases/tags/{tag}"


def _headers(token):
    headers = {
        "Accept": "application/vnd.github+json",
        "User-Agent": "PySME-release-asset-downloader",
    }
    if token:
        headers["Authorization"] = f"Bearer {token}"
    return headers


def fetch_release(repo, tag, token):
    url = API_TEMPLATE.format(repo=repo, tag=tag)
    req = Request(url, headers=_headers(token))
    with urlopen(req) as resp:
        return json.loads(resp.read().decode("utf-8"))


def download_file(url, destination, token):
    req = Request(url, headers=_headers(token))
    with urlopen(req) as resp, open(destination, "wb") as f:
        while True:
            chunk = resp.read(1024 * 1024)
            if not chunk:
                break
            f.write(chunk)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", default="SpectroscopyMadeEasy/PySME", help="GitHub repository (owner/name)")
    parser.add_argument("--tag", required=True, help="Release tag, e.g. v1.0.0")
    parser.add_argument("--output", default=None, help="Output directory. Default: ./release_assets/<tag>")
    parser.add_argument("--force", action="store_true", help="Re-download files even if already present")
    parser.add_argument("--token-env", default="GITHUB_TOKEN", help="Environment variable for optional GitHub token")
    return parser.parse_args()


def main():
    args = parse_args()
    token = os.environ.get(args.token_env)

    outdir = Path(args.output) if args.output else Path("release_assets") / args.tag
    outdir.mkdir(parents=True, exist_ok=True)

    try:
        release = fetch_release(args.repo, args.tag, token)
    except HTTPError as exc:
        print(f"[error] GitHub API returned HTTP {exc.code}: {exc.reason}", file=sys.stderr)
        return 1
    except URLError as exc:
        print(f"[error] Could not reach GitHub API: {exc.reason}", file=sys.stderr)
        return 1

    assets = release.get("assets", [])
    if not assets:
        print("[warn] No assets found on this release.")

    manifest = {
        "repo": args.repo,
        "tag": args.tag,
        "release_name": release.get("name"),
        "published_at": release.get("published_at"),
        "html_url": release.get("html_url"),
        "assets": [],
    }

    for asset in assets:
        name = asset["name"]
        url = asset["browser_download_url"]
        destination = outdir / name

        if destination.exists() and not args.force:
            print(f"[skip] {name} (already exists)")
        else:
            print(f"[download] {name}")
            try:
                download_file(url, destination, token)
            except HTTPError as exc:
                print(f"[error] Failed to download {name}: HTTP {exc.code} {exc.reason}", file=sys.stderr)
                return 1
            except URLError as exc:
                print(f"[error] Failed to download {name}: {exc.reason}", file=sys.stderr)
                return 1

        manifest["assets"].append(
            {
                "name": name,
                "size": asset.get("size"),
                "digest": asset.get("digest"),
                "download_url": url,
                "local_path": str(destination),
            }
        )

    manifest_path = outdir / "manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n")
    print(f"[done] Downloaded {len(assets)} assets to {outdir}")
    print(f"[done] Wrote manifest: {manifest_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
