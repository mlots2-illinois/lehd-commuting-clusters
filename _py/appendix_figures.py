"""
appendix_figures.py -- schematic figures for the appendix: the snake-scan window diagram, the Leiden resolution calibration figure, and the administrative-units-versus-places schematic.

Called from gen/fig_snakescan.do and gen/fig_admin_places.do as:
    python appendix_figures.py $tables/leiden_resolution_sweep_tract.csv <figures_dir> snakescan|admin
    also accepted: all | calibration as the third argument, and --picked-gamma G for the calibration figure
Reads:  <sweep_csv> (calibration mode only; the positional is required in every mode so the Stata call stays fixed)
Writes: <figures_dir>/appendix_snakescan.png; <figures_dir>/appendix_admin_vs_places.png; <figures_dir>/appendix_calibration.png (calibration mode)
Notes:  the calibration mode re-derives the selected gamma from the sweep (the largest gamma within 0.05 of peak stability among non-degenerate resolutions) and warns when --picked-gamma disagrees.
        The snakescan and admin figures are drawn from constants, not from data.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

C_BLUE = "#3b5bdb"
C_GREEN = "#2f9e44"
C_MAROON = "#c92a2a"
C_GRAY = "#6b6b6b"
C_BG = "#f0f0f0"

plt.rcParams.update({
    "font.size": 11,
    "font.family": "serif",
    "axes.linewidth": 0.8,
})

def fig_snakescan(out_path: Path) -> None:
    fig, (ax1, ax2) = plt.subplots(
        2, 1, figsize=(8.2, 5.4), gridspec_kw={"height_ratios": [1.25, 1.0]}
    )

    nrows, ncols = 4, 8
    xs, ys, order = [], [], []
    k = 0
    for r in range(nrows):
        cols = range(ncols) if r % 2 == 0 else range(ncols - 1, -1, -1)
        for c in cols:
            xs.append(c)
            ys.append(nrows - 1 - r)
            order.append(k)
            k += 1
    xs, ys = np.array(xs, float), np.array(ys, float)

    ax1.scatter(xs, ys, s=170, facecolor=C_BG, edgecolor="black",
                zorder=3, linewidths=0.9)
    for i in range(len(xs) - 1):
        ax1.annotate(
            "", xy=(xs[i + 1], ys[i + 1]), xytext=(xs[i], ys[i]),
            arrowprops=dict(arrowstyle="-|>", color=C_GRAY, lw=1.4,
                            shrinkA=9, shrinkB=9),
            zorder=2,
        )
    ax1.text((ncols - 1) / 2, nrows - 0.25,
             "each square = one tract; arrows show the snake-scan ordering",
             ha="center", va="bottom", fontsize=9.5, color=C_GRAY)
    ax1.text(-0.85, nrows - 1.5, "latitude\nbands", ha="center", va="center",
             fontsize=9.5, color=C_GRAY, rotation=90)
    ax1.set_xlim(-1.4, ncols - 0.3)
    ax1.set_ylim(-0.7, nrows + 0.2)
    ax1.set_title("(a)  A two-dimensional field is folded into a one-dimensional order",
                  fontsize=11, loc="left")
    ax1.axis("off")

    n = 16
    px = np.arange(n)
    ax2.scatter(px, np.zeros(n), s=150, facecolor=C_BG, edgecolor="black",
                zorder=3, linewidths=0.9)
    ax2.plot([0, n - 1], [0, 0], color=C_GRAY, lw=1.0, zorder=1)

    def window(lo, hi, color, y, label, labtop=True):
        box = FancyBboxPatch(
            (lo - 0.45, y - 0.42), (hi - lo) + 0.9, 0.84,
            boxstyle="round,pad=0.02,rounding_size=0.18",
            linewidth=1.8, edgecolor=color, facecolor=color, alpha=0.16,
            zorder=2,
        )
        ax2.add_patch(box)
        ly = y + (0.62 if labtop else -0.62)
        ax2.text((lo + hi) / 2, ly, label, ha="center",
                 va="bottom" if labtop else "top", color=color, fontsize=10)

    window(0, 8, C_BLUE, 0.0, "window 1", labtop=True)
    window(6, 14, C_GREEN, 0.0, "window 2", labtop=False)

    for u in (6, 7, 8):
        ax2.scatter([u], [0], s=160, facecolor=C_MAROON, edgecolor="black",
                    zorder=4, linewidths=0.9, alpha=0.85)
    ax2.annotate(
        "overlap: clustered in BOTH windows\n(used to stitch the seam)",
        xy=(7, 0), xytext=(7, 1.55), ha="center", va="bottom",
        fontsize=9.5, color=C_MAROON,
        arrowprops=dict(arrowstyle="-|>", color=C_MAROON, lw=1.3),
    )
    ax2.set_xlim(-1.2, n - 0.2)
    ax2.set_ylim(-1.5, 2.4)
    ax2.set_title("(b)  Fixed-width windows slide along the order with a shorter stride, so they overlap",
                  fontsize=11, loc="left")
    ax2.axis("off")

    fig.tight_layout(h_pad=1.6)
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out_path}")

def fig_calibration(sweep_csv: Path, out_path: Path,
                    picked_gamma: float | None = None) -> None:
    df = pd.read_csv(sweep_csv).sort_values("resolution").reset_index(drop=True)
    g = df["resolution"].to_numpy(float)
    stab = df["mean_stability"].to_numpy(float)
    nclu = df["n_clusters_est_national"].to_numpy(float)
    sing = df["mean_singleton_frac"].to_numpy(float)
    lfrac = df["mean_largest_frac"].to_numpy(float)

    SING_MAX, DEGEN, TOL = 0.50, 0.50, 0.05
    ok = (sing < SING_MAX) & (lfrac < DEGEN)
    peak = stab[ok].max()
    near = ok & (stab >= peak - TOL)
    derived_g = g[near].max()
    if picked_gamma is None:
        pick_g = derived_g
    else:
        pick_g = picked_gamma
        if abs(derived_g - picked_gamma) > 1e-9:
            print(f"WARNING: --picked-gamma {picked_gamma:g} disagrees with the "
                  f"re-derived selection {derived_g:g} from {sweep_csv.name}; "
                  f"drawing the supplied value; check that the sweep table and "
                  f"the persisted pick come from the same run.", file=sys.stderr)

    fig, axL = plt.subplots(figsize=(8.2, 4.6))
    axR = axL.twinx()

    axL.axhspan(peak - TOL, 1.001, color=C_BLUE, alpha=0.07, zorder=0)
    axL.text(g.min(), peak - TOL, "  reproducible plateau (within tol. of peak)",
             va="top", ha="left", fontsize=8.5, color=C_BLUE)

    lS, = axL.plot(g, stab, "-o", color=C_BLUE, lw=1.8, ms=5,
                   label="cross-seed reproducibility")
    axL.set_ylabel("cross-seed reproducibility\n(mean pairwise agreement)",
                   color=C_BLUE)
    axL.tick_params(axis="y", labelcolor=C_BLUE)
    span = max(stab.max() - stab.min(), 1e-3)
    y0 = stab.min() - 0.15 * span
    y1 = max(stab.max(), 1.0) + 0.08 * span
    axL.set_ylim(y0, y1)

    lN, = axR.plot(g, nclu, "--s", color=C_GREEN, lw=1.6, ms=4,
                   label="estimated cluster count")
    axR.set_yscale("log")
    axR.set_ylabel("estimated cluster count (log)", color=C_GREEN)
    axR.tick_params(axis="y", labelcolor=C_GREEN)

    axL.set_xscale("log")
    axL.set_xlabel(r"Leiden resolution  $\gamma$  (log scale)")

    axL.axvline(pick_g, color="black", lw=1.4, ls=":")
    axL.annotate(
        rf"selected $\gamma = {pick_g:g}$" + "\n(finest reproducible scale)",
        xy=(pick_g, y0 + 0.30 * (y1 - y0)), xytext=(pick_g * 0.18, y0 + 0.22 * (y1 - y0)),
        fontsize=9.5, ha="left",
        arrowprops=dict(arrowstyle="-|>", color="black", lw=1.2),
    )

    deg = ~ok
    if deg.any():
        gd, sd = g[deg], stab[deg]
        axL.scatter(gd, sd, marker="x", s=90, color=C_MAROON, zorder=5,
                    linewidths=2.0)
        gi = int(np.where(deg)[0][-1])
        axL.annotate(
            "rejected:\nall singletons\n(trivially reproducible)",
            xy=(g[gi], stab[gi]), xytext=(g[gi] * 0.95, y0 + 0.60 * (y1 - y0)),
            fontsize=8.8, ha="right", color=C_MAROON,
            arrowprops=dict(arrowstyle="-|>", color=C_MAROON, lw=1.1),
        )

    axL.legend(handles=[lS, lN], loc="lower left", fontsize=9, framealpha=0.9)
    axL.set_title("Selecting the Leiden resolution by reproducibility, not by matching k",
                  fontsize=11, loc="left")
    fig.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out_path} (picked gamma={pick_g:g}, peak stab={peak:.3f})")

def fig_admin_vs_places(out_path: Path) -> None:
    region = np.array([
        [0, 0, 0, 1, 1, 1],
        [0, 0, 0, 1, 1, 1],
        [0, 0, 1, 1, 1, 1],
        [2, 2, 2, 2, 1, 1],
        [2, 2, 2, 2, 2, 1],
        [2, 2, 2, 2, 2, 2],
    ])
    cols = {0: C_BLUE, 1: C_GREEN, 2: C_MAROON}
    nr, nc = region.shape

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8.2, 4.3))

    def draw_grid(ax, colored):
        for r in range(nr):
            for c in range(nc):
                y = nr - 1 - r
                fc = (cols[region[r, c]] if colored else C_BG)
                alpha = 0.5 if colored else 1.0
                ax.add_patch(mpatches.Rectangle(
                    (c, y), 1, 1, facecolor=fc, alpha=alpha,
                    edgecolor="white" if colored else C_GRAY,
                    linewidth=1.0 if colored else 0.7))
        ax.set_xlim(-0.1, nc + 0.1)
        ax.set_ylim(-0.1, nr + 0.1)
        ax.set_aspect("equal")
        ax.axis("off")

    draw_grid(ax1, colored=False)
    ax1.plot([nc / 2, nc / 2], [0, nr], color="black", lw=2.4, ls="--")
    ax1.plot([0, nc], [nr / 2, nr / 2], color="black", lw=2.4, ls="--")
    ax1.set_title("(a)  Administrative units (input):\nboundaries drawn for record-keeping",
                  fontsize=10.5, loc="left")

    draw_grid(ax2, colored=True)
    ax2.set_title("(b)  Delineated places (output):\nboundaries that follow commuting",
                  fontsize=10.5, loc="left")

    fig.tight_layout(w_pad=2.0)
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out_path}")

def main() -> None:
    parser = argparse.ArgumentParser(
        prog="appendix_figures.py",
        description="Render the Appendix A schematic figures.",
        usage=("appendix_figures.py <sweep_csv> <out_dir> "
               "[all|snakescan|calibration|admin] [--picked-gamma G]\n"
               "  sweep_csv is required even for modes that ignore it "
               "(snakescan/admin); kept positional for Stata-call stability."))
    parser.add_argument("sweep_csv",
                        help="leiden_resolution_sweep_<level>.csv (used by the "
                             "calibration figure; still required for other modes)")
    parser.add_argument("out_dir", help="output directory for the PNGs")
    parser.add_argument("which", nargs="?", default="all",
                        choices=["all", "snakescan", "calibration", "admin"],
                        help="which figure(s) to render (default: all)")
    parser.add_argument("--picked-gamma", type=float, default=None,
                        help="canonical selected gamma to mark in the "
                             "calibration figure (the internal re-derivation "
                             "becomes a cross-check)")
    args = parser.parse_args()

    sweep_csv = Path(args.sweep_csv)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    if args.which in ("all", "snakescan"):
        fig_snakescan(out_dir / "appendix_snakescan.png")
    if args.which in ("all", "calibration"):
        fig_calibration(sweep_csv, out_dir / "appendix_calibration.png",
                        picked_gamma=args.picked_gamma)
    if args.which in ("all", "admin"):
        fig_admin_vs_places(out_dir / "appendix_admin_vs_places.png")

if __name__ == "__main__":
    main()
