#!/usr/bin/env python3
"""Generate the two PySME v1.2 documentation figures from audit artifacts."""

from __future__ import annotations

import json
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import FuncFormatter


ROOT = Path(__file__).resolve().parents[1]
ANALYSIS = ROOT / "analysis"
OUTPUT = ROOT / "docs" / "_static" / "v120"

HISTORICAL_COLOR = "#6B7280"
V120_COLOR = "#2F6F9F"
RESIDUAL_COLOR = "#B45309"
GRID_COLOR = "#D1D5DB"
TEXT_COLOR = "#1F2937"


def configure_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 10,
            "axes.labelsize": 10,
            "axes.titlesize": 11,
            "axes.titleweight": "semibold",
            "axes.edgecolor": "#9CA3AF",
            "axes.labelcolor": TEXT_COLOR,
            "xtick.color": TEXT_COLOR,
            "ytick.color": TEXT_COLOR,
            "text.color": TEXT_COLOR,
            "legend.frameon": False,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "savefig.facecolor": "white",
            "savefig.bbox": "tight",
        }
    )


def load_release_spectrum() -> dict[str, object]:
    benchmark_path = ANALYSIS / "v120_release_benchmark.json"
    benchmark = json.loads(benchmark_path.read_text())
    case = benchmark["cases"]["cool_metal_rich_dwarf_10A"]
    historical = case["profiles"]["H"]
    current = case["profiles"]["C1B1"]

    old_array = ROOT / historical["runs"][0]["array_file"]
    new_array = ROOT / current["runs"][0]["array_file"]
    with np.load(old_array) as old_data, np.load(new_array) as new_data:
        wave = old_data["wave"].copy()
        old_flux = old_data["flux"].copy()
        new_wave = new_data["wave"].copy()
        new_flux = new_data["flux"].copy()

    if not np.array_equal(wave, new_wave):
        raise RuntimeError("Historical and v1.2 output grids differ")

    residual = new_flux - old_flux
    old_time = float(historical["median_sec"])
    new_time = float(current["median_sec"])
    return {
        "wave": wave,
        "old_flux": old_flux,
        "new_flux": new_flux,
        "residual": residual,
        "old_time": old_time,
        "new_time": new_time,
        "speedup": old_time / new_time,
        "max_abs": float(np.max(np.abs(residual))),
        "rms": float(np.sqrt(np.mean(residual**2))),
        "benchmark_path": benchmark_path,
        "old_array": old_array,
        "new_array": new_array,
    }


def load_candidate_work() -> dict[str, object]:
    benchmark_path = ANALYSIS / "batched_rkints_benchmark.json"
    benchmark = json.loads(benchmark_path.read_text())
    model = benchmark["models"]["cool_metal_rich_dwarf"]
    legacy_nodes = int(model["legacy"]["transfer_points"])
    current_nodes = int(model["batched"]["transfer_points"])

    timing_path = ANALYSIS / "batched_rkints_cool_timing.log"
    pattern = re.compile(
        r"intervals\s+: calls=(\d+), candidates=(\d+), "
        r"full_scan=(\d+), reduction=([0-9.]+)%"
    )
    generations: list[tuple[int, int, int, float]] = []
    for match in pattern.finditer(timing_path.read_text()):
        calls, candidates, full_scan, reduction = match.groups()
        calls_i = int(calls)
        if calls_i in {5025, 5024}:
            generations.append(
                (calls_i, int(candidates), int(full_scan), float(reduction))
            )
    if len(generations) != 2:
        raise RuntimeError("Could not identify both cool-star batched generations")

    current_visits = sum(item[1] for item in generations)
    legacy_visits = sum(item[2] for item in generations)
    reduction = 1.0 - current_visits / legacy_visits
    return {
        "legacy_nodes": legacy_nodes,
        "current_nodes": current_nodes,
        "legacy_visits": legacy_visits,
        "current_visits": current_visits,
        "reduction": reduction,
        "benchmark_path": benchmark_path,
        "timing_path": timing_path,
    }


def add_panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(
        -0.10,
        1.04,
        label,
        transform=ax.transAxes,
        fontsize=12,
        fontweight="semibold",
        va="bottom",
    )


def save_figure(fig: plt.Figure, stem: str) -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT / f"{stem}.png", dpi=200)
    fig.savefig(OUTPUT / f"{stem}.pdf")
    plt.close(fig)


def make_spectral_comparison(data: dict[str, object]) -> None:
    wave = np.asarray(data["wave"])
    old_flux = np.asarray(data["old_flux"])
    new_flux = np.asarray(data["new_flux"])
    residual = np.asarray(data["residual"])

    fig, (ax_flux, ax_resid) = plt.subplots(
        2,
        1,
        figsize=(8.0, 6.0),
        sharex=True,
        gridspec_kw={"height_ratios": [2.2, 1.0], "hspace": 0.08},
    )

    ax_flux.plot(
        wave,
        old_flux,
        color=HISTORICAL_COLOR,
        linewidth=1.6,
        linestyle=(0, (4, 2)),
        label="Historical",
    )
    ax_flux.plot(
        wave,
        new_flux,
        color=V120_COLOR,
        linewidth=1.25,
        label="PySME v1.2",
    )
    ax_flux.set_ylabel("Normalized flux")
    ax_flux.set_title("Cool metal-rich dwarf, 5195–5205 Å", loc="left")
    ax_flux.legend(loc="lower right", ncol=2)
    ax_flux.grid(axis="y", color=GRID_COLOR, linewidth=0.6, alpha=0.65)
    ax_flux.text(
        0.985,
        0.965,
        "\n".join(
            [
                f"Historical: {data['old_time']:.1f} s",
                f"PySME v1.2: {data['new_time']:.3f} s",
                f"Speed-up: {data['speedup']:.1f}×",
            ]
        ),
        transform=ax_flux.transAxes,
        ha="right",
        va="top",
        bbox={
            "boxstyle": "round,pad=0.35",
            "facecolor": "white",
            "edgecolor": "#D1D5DB",
            "alpha": 0.94,
        },
    )
    add_panel_label(ax_flux, "A")

    resid_limit = max(1.55e-9, 1.08 * float(np.max(np.abs(residual))))
    ax_resid.axhline(0.0, color="#9CA3AF", linewidth=0.8)
    ax_resid.plot(wave, residual, color=RESIDUAL_COLOR, linewidth=1.2)
    ax_resid.set_ylim(-resid_limit, resid_limit)
    ax_resid.set_xlabel("Wavelength (Å)")
    ax_resid.set_ylabel(r"$\Delta F$")
    ax_resid.ticklabel_format(axis="y", style="sci", scilimits=(-2, 2))
    ax_resid.grid(axis="y", color=GRID_COLOR, linewidth=0.6, alpha=0.65)
    ax_resid.text(
        0.985,
        0.92,
        "\n".join(
            [
                rf"max $|\Delta F|$ = {compact_scientific(data['max_abs'])}",
                rf"RMS = {compact_scientific(data['rms'])}",
            ]
        ),
        transform=ax_resid.transAxes,
        ha="right",
        va="top",
        bbox={
            "boxstyle": "round,pad=0.30",
            "facecolor": "white",
            "edgecolor": "#D1D5DB",
            "alpha": 0.94,
        },
    )
    add_panel_label(ax_resid, "B")

    for ax in (ax_flux, ax_resid):
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.margins(x=0)

    save_figure(fig, "v120_spectral_comparison")


def compact_scientific(value: float) -> str:
    exponent = int(np.floor(np.log10(value)))
    coefficient = value / 10**exponent
    return rf"${coefficient:.2f}\times10^{{{exponent}}}$"


def make_candidate_work(data: dict[str, object]) -> None:
    labels = ["Historical", "PySME v1.2"]
    colors = [HISTORICAL_COLOR, V120_COLOR]
    fig, (ax_nodes, ax_visits) = plt.subplots(1, 2, figsize=(9.0, 4.3))

    nodes = [int(data["legacy_nodes"]), int(data["current_nodes"])]
    bars_nodes = ax_nodes.bar(labels, nodes, color=colors, width=0.58)
    ax_nodes.set_title("Physical transfer wavelengths")
    ax_nodes.set_ylabel("Transfer nodes")
    ax_nodes.set_ylim(0, max(nodes) * 1.22)
    ax_nodes.yaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{x:,.0f}"))
    ax_nodes.grid(axis="y", color=GRID_COLOR, linewidth=0.6, alpha=0.65)
    for bar, value in zip(bars_nodes, nodes):
        ax_nodes.text(
            bar.get_x() + bar.get_width() / 2,
            value + max(nodes) * 0.035,
            f"{value:,}",
            ha="center",
            va="bottom",
            fontweight="semibold",
        )
    add_panel_label(ax_nodes, "A")

    visits = [int(data["legacy_visits"]), int(data["current_visits"])]
    bars_visits = ax_visits.bar(labels, visits, color=colors, width=0.58)
    ax_visits.set_title("Candidate-line visits")
    ax_visits.set_ylabel("Candidate visits (log scale)")
    ax_visits.set_yscale("log")
    ax_visits.set_ylim(1e6, 1.1e10)
    ax_visits.grid(axis="y", which="major", color=GRID_COLOR, linewidth=0.6, alpha=0.65)
    for bar, value in zip(bars_visits, visits):
        ax_visits.text(
            bar.get_x() + bar.get_width() / 2,
            value * 1.20,
            compact_scientific(value),
            ha="center",
            va="bottom",
            fontweight="semibold",
        )
    ax_visits.text(
        0.97,
        0.94,
        f"{100 * data['reduction']:.2f}% reduction",
        transform=ax_visits.transAxes,
        ha="right",
        va="top",
        bbox={
            "boxstyle": "round,pad=0.30",
            "facecolor": "white",
            "edgecolor": "#D1D5DB",
            "alpha": 0.94,
        },
    )
    add_panel_label(ax_visits, "B")

    for ax in (ax_nodes, ax_visits):
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.tick_params(axis="x", length=0)

    fig.subplots_adjust(wspace=0.34)
    save_figure(fig, "v120_candidate_work_reduction")


def main() -> None:
    configure_style()
    spectral = load_release_spectrum()
    candidate = load_candidate_work()
    make_spectral_comparison(spectral)
    make_candidate_work(candidate)
    print(
        json.dumps(
            {
                "spectral": {
                    "old_time": spectral["old_time"],
                    "new_time": spectral["new_time"],
                    "speedup": spectral["speedup"],
                    "max_abs": spectral["max_abs"],
                    "rms": spectral["rms"],
                },
                "candidate_work": {
                    "legacy_nodes": candidate["legacy_nodes"],
                    "current_nodes": candidate["current_nodes"],
                    "legacy_visits": candidate["legacy_visits"],
                    "current_visits": candidate["current_visits"],
                    "reduction": candidate["reduction"],
                },
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
