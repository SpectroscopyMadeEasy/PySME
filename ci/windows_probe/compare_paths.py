#!/usr/bin/env python3
"""Compare Unicode-path synthesis results with the ASCII-path baseline."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

from annotations import emit


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--ascii", type=Path, required=True)
    parser.add_argument("--unicode", type=Path, required=True)
    args = parser.parse_args()
    failed = False
    with np.load(args.ascii, allow_pickle=False) as baseline, np.load(args.unicode, allow_pickle=False) as candidate:
        for name in ("solar_fixed", "hydrogen_fixed"):
            wave_a = baseline[f"{name}__wave"]
            wave_b = candidate[f"{name}__wave"]
            flux_a = baseline[f"{name}__flux"]
            flux_b = candidate[f"{name}__flux"]
            if not np.array_equal(wave_a, wave_b):
                emit("error", f"Unicode path {name}", "wavelength grids differ")
                failed = True
                continue
            diff = np.asarray(flux_b - flux_a, dtype=float)
            max_abs = float(np.max(np.abs(diff)))
            rms = float(np.sqrt(np.mean(diff**2)))
            emit("notice", f"Unicode path {name}", f"pixels={diff.size} max={max_abs:.9g} rms={rms:.9g}")
            if max_abs > 1e-12:
                emit("error", f"Unicode path mismatch {name}", f"max={max_abs:.9g} > 1e-12")
                failed = True
    return int(failed)


if __name__ == "__main__":
    raise SystemExit(main())
