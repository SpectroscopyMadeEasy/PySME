#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from __future__ import annotations

import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import numpy as np
from scipy.ndimage import binary_fill_holes

from pysme.atmosphere.savfile import SavFile
from pysme.large_file_storage import setup_lfs


ATMOSPHERE_FILES = [
    "marcs2012.sav",
    "marcs2012p_t0.0.sav",
    "marcs2012p_t1.0.sav",
    "marcs2012p_t2.0.sav",
    "marcs2012s_t1.0.sav",
    "marcs2012s_t2.0.sav",
    "marcs2012s_t5.0.sav",
    "marcs2012t00cooldwarfs.sav",
    "marcs2012t01cooldwarfs.sav",
    "marcs2012t02cooldwarfs.sav",
    "atlas12.sav",
    "atlas9_vmic0.0.sav",
    "atlas9_vmic2.0.sav",
    "ll_vmic2.0.sav",
]


def build_filled_mask(grid):
    teff_values = np.sort(np.unique(np.asarray(grid.teff, dtype=float)))
    logg_values = np.sort(np.unique(np.asarray(grid.logg, dtype=float)))
    teff_index = {value: idx for idx, value in enumerate(teff_values)}
    logg_index = {value: idx for idx, value in enumerate(logg_values)}
    monh_values = np.sort(np.unique(np.asarray(grid.monh, dtype=float)))[::-1]

    masks = {}
    for monh in monh_values:
        select = np.isclose(grid.monh, monh)
        occupied = np.zeros((len(logg_values), len(teff_values)), dtype=bool)
        for teff, logg in zip(
            np.asarray(grid.teff[select], dtype=float),
            np.asarray(grid.logg[select], dtype=float),
        ):
            occupied[logg_index[logg], teff_index[teff]] = True

        row_fill = occupied.copy()
        for i in range(row_fill.shape[0]):
            cols = np.flatnonzero(row_fill[i])
            if cols.size > 0:
                row_fill[i, cols.min() : cols.max() + 1] = True

        col_fill = row_fill.copy()
        for j in range(col_fill.shape[1]):
            rows = np.flatnonzero(col_fill[:, j])
            if rows.size > 0:
                col_fill[rows.min() : rows.max() + 1, j] = True

        masks[monh] = binary_fill_holes(col_fill)

    return monh_values, teff_values, logg_values, masks


def boundary_segments(mask, teff_edges, logg_edges):
    segments = []
    nrows, ncols = mask.shape
    for i in range(nrows):
        for j in range(ncols):
            if not mask[i, j]:
                continue
            x0, x1 = teff_edges[j], teff_edges[j + 1]
            y0, y1 = logg_edges[i], logg_edges[i + 1]
            if i == 0 or not mask[i - 1, j]:
                segments.append([(x0, y0), (x1, y0)])
            if i == nrows - 1 or not mask[i + 1, j]:
                segments.append([(x0, y1), (x1, y1)])
            if j == 0 or not mask[i, j - 1]:
                segments.append([(x0, y0), (x0, y1)])
            if j == ncols - 1 or not mask[i, j + 1]:
                segments.append([(x1, y0), (x1, y1)])
    return segments


def plot_grid(source, output):
    _, lfs_atmo, _ = setup_lfs()
    grid_path = lfs_atmo.get(source)
    grid = SavFile(grid_path, source=source, lfs=None)

    monh_values, teff_values, logg_values, masks = build_filled_mask(grid)
    teff_step = float(np.min(np.diff(teff_values))) if len(teff_values) > 1 else 100.0
    logg_step = float(np.min(np.diff(logg_values))) if len(logg_values) > 1 else 0.5
    teff_edges = np.concatenate(
        (
            [teff_values[0] - teff_step / 2],
            (teff_values[:-1] + teff_values[1:]) / 2,
            [teff_values[-1] + teff_step / 2],
        )
    )
    logg_edges = np.concatenate(
        (
            [logg_values[0] - logg_step / 2],
            (logg_values[:-1] + logg_values[1:]) / 2,
            [logg_values[-1] + logg_step / 2],
        )
    )

    ncols = 4
    nrows = math.ceil(len(monh_values) / ncols)
    fig, axes = plt.subplots(
        nrows,
        ncols,
        figsize=(4 * ncols, 3.1 * nrows),
        squeeze=False,
        sharex=True,
        sharey=True,
    )

    teff_min = float(np.nanmin(grid.teff))
    teff_max = float(np.nanmax(grid.teff))
    logg_edge_min = float(logg_edges[0])
    logg_edge_max = float(logg_edges[-1])
    add_boundary = source == "marcs2012.sav"

    for ax, monh in zip(axes.ravel(), monh_values):
        select = np.isclose(grid.monh, monh)
        teff = np.asarray(grid.teff[select], dtype=float)
        logg = np.asarray(grid.logg[select], dtype=float)
        ax.scatter(teff, logg, s=8, alpha=0.8, color="#2c7fb8", zorder=3)
        if add_boundary:
            segments = boundary_segments(masks[monh], teff_edges, logg_edges)
            ax.add_collection(
                LineCollection(segments, colors="#d7301f", linewidths=2.0, zorder=4)
            )
        ax.set_title(f"[M/H] = {monh:g}")
        ax.set_xlim(teff_max + teff_step, teff_min - teff_step)
        ax.set_ylim(logg_edge_max + 0.1, logg_edge_min - 0.1)
        ax.grid(alpha=0.2)

    for ax in axes.ravel()[len(monh_values) :]:
        ax.axis("off")
    for ax in axes[-1, :]:
        ax.set_xlabel("Teff [K]")
    for ax in axes[:, 0]:
        ax.set_ylabel("log g")

    fig.tight_layout()
    fig.savefig(output, dpi=180)
    plt.close(fig)


def main():
    docs_dir = Path(__file__).resolve().parents[1] / "docs" / "img" / "atmosphere"
    docs_dir.mkdir(parents=True, exist_ok=True)
    for source in ATMOSPHERE_FILES:
        output = docs_dir / f"{source.replace('.sav', '')}_grid.png"
        print(f"Generating {output.name}")
        plot_grid(source, output)


if __name__ == "__main__":
    main()
