#!/usr/bin/env python3
"""Plot the principal residuals from the RKINTS second-pruning audit."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent
DATA = np.load(ROOT / "rkints_second_pruning_audit.npz")


def central(wave, lo, hi):
    return (wave >= lo) & (wave <= hi)


fig, axes = plt.subplots(3, 1, figsize=(10, 8.2), constrained_layout=True)

# The fixed fine grid isolates active-mask semantics from adaptive-grid errors.
w = DATA["metal_poor_dwarf_fine_wave"]
m = central(w, 5195.0, 5205.0)
a = DATA["metal_poor_dwarf_A_almax_only_fine_flux"]
l = DATA["metal_poor_dwarf_L_final_mask_fine_flux"]
f = DATA["metal_poor_dwarf_F_1e-6_fine_flux"]
axes[0].plot(w[m], 1e3 * (l[m] - a[m]), lw=0.9, color="tab:red")
axes[0].axhline(0, color="0.4", lw=0.6)
axes[0].set_ylabel(r"$(F_L-F_A)\times10^3$")
axes[0].set_title("Metal-poor dwarf: isolated effect of legacy MARK=2 pruning")
axes[0].annotate(
    "Ti I 5201.0814 Å",
    xy=(5201.0814, np.interp(5201.0814, w[m], 1e3 * (l[m] - a[m]))),
    xytext=(5199.7, 1.45),
    arrowprops={"arrowstyle": "->", "lw": 0.8},
    fontsize=9,
)

axes[1].plot(w[m], 1e3 * (l[m] - f[m]), lw=0.9, label="Legacy final mask − F", color="tab:red")
axes[1].plot(w[m], 1e3 * (a[m] - f[m]), lw=0.9, label="ALMAX-only − F", color="tab:blue")
axes[1].axhline(0, color="0.4", lw=0.6)
axes[1].set_ylabel(r"residual $\times10^3$")
axes[1].set_title("Same fine grid against minimally pruned F (range floor $10^{-6}$)")
axes[1].legend(loc="upper right", frameon=False, fontsize=9)

# Only the weak-line decision order is reversed; opacity accumulation order is fixed.
wc = DATA["extra_metal_poor_dwarf_clean_atomic_wave"]
mc = central(wc, 5300.0, 5310.0)
forward = DATA["extra_metal_poor_dwarf_clean_atomic_audit-forward_flux"]
reverse = DATA["extra_metal_poor_dwarf_clean_atomic_audit-reverse_flux"]
axes[2].plot(wc[mc], 1e3 * (reverse[mc] - forward[mc]), lw=0.9, color="tab:purple")
axes[2].axhline(0, color="0.4", lw=0.6)
axes[2].set_ylabel(r"$(F_{rev}-F_{fwd})\times10^3$")
axes[2].set_xlabel("Wavelength (Å)")
axes[2].set_title("Clean atomic window: decision-order dependence (169 vs 171 pruned lines)")

for ax in axes:
    ax.grid(alpha=0.18, linewidth=0.5)
    ax.margins(x=0)

out = ROOT / "rkints_second_pruning_residuals.png"
fig.savefig(out, dpi=180)
print(out)
