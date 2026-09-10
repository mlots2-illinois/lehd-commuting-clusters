"""
calib_dispersion_lines.py -- three-panel line figure of D dispersion (mean, SD, share above 0.95) against the distance cap, one line per alpha, per level.

Called from gen/fig_calib_dispersion.do as:
    python calib_dispersion_lines.py <tables_dir> <figures_dir>
    also accepted: one or more levels after <figures_dir>; without them both tract and county are drawn
Reads:  <tables_dir>/calib_dispersion_<level>.csv
Writes: <figures_dir>/calib_dispersion_lines_<level>.png
Notes:  at most seven alpha values are drawn, evenly spaced through the grid, so the lines stay distinguishable. The cap axis is logarithmic.
"""

import os
import sys
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

N_LINES = 7

def alpha_subset(alphas):
    alphas = sorted(set(round(a, 6) for a in alphas))
    if len(alphas) <= N_LINES:
        return alphas
    idx = [round(i * (len(alphas) - 1) / (N_LINES - 1)) for i in range(N_LINES)]
    return sorted({alphas[i] for i in idx})

def plot_level(tables_dir, figures_dir, level):
    csv = f"{tables_dir}/calib_dispersion_{level}.csv"
    if not os.path.exists(csv):
        print(f"  (skip {level}: {csv} not found)")
        return
    df = pd.read_csv(csv)
    lines = alpha_subset(df["alpha"].unique())
    cmap = plt.cm.viridis
    panels = [
        ("d_mean", "Mean $D_{ij}$", "interior = informative; near 0/1 = collapsed"),
        ("d_sd", "SD of $D_{ij}$", "want large = informative spread"),
        ("d_frac_hi", "Share $D_{ij} > 0.95$", "want low = not collapsed to 1"),
    ]
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.6))
    for ax, (col, ylab, sub) in zip(axes, panels):
        for i, a in enumerate(lines):
            g = df[abs(df["alpha"] - a) < 1e-6].sort_values("cap_km")
            if g.empty:
                continue
            ax.plot(g["cap_km"], g[col], marker="o", ms=3, lw=1.5,
                    color=cmap(i / max(len(lines) - 1, 1)),
                    label=f"α={a:g}")
        ax.set_xscale("log")
        ax.set_xlabel("distance cap  $\\bar d$  (km)")
        ax.set_ylabel(ylab)
        ax.set_title(ylab, fontsize=10)
        ax.grid(alpha=0.25)
    axes[0].legend(fontsize=8, loc="best", title="flow weight α")
    fig.suptitle(f"Dispersion of $D_{{ij}}$ over the (α × cap) grid — {level}",
                 fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    out = f"{figures_dir}/calib_dispersion_lines_{level}.png"
    fig.savefig(out, dpi=140)
    plt.close(fig)
    print("  ->", out)

if __name__ == "__main__":
    if len(sys.argv) < 3:
        sys.exit("usage: calib_dispersion_lines.py <tables_dir> <figures_dir> "
                 "[level ...]")
    tables_dir, figures_dir = sys.argv[1], sys.argv[2]
    levels = sys.argv[3:] or ["tract", "county"]
    for lvl in levels:
        plot_level(tables_dir, figures_dir, lvl)
