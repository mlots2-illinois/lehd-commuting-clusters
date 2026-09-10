"""
leiden_resolution_lines.py -- figures of the Leiden resolution sweep: a two-panel compact view (cluster count; largest-community and singleton fractions) or a four-panel diagnostic view.

Called from calib/calib_02_resolution.do as:
    python leiden_resolution_lines.py <tables_dir> <figures_dir> <level> --picked-gamma <res> --compact
and from gen/fig_gamma_fine.do as:
    python leiden_resolution_lines.py <tables_dir> <figures_dir> <level> [--picked-gamma <res>]
Reads:  <tables_dir>/leiden_resolution_sweep_<level>.csv
Writes: <figures_dir>/leiden_resolution_lines_<level>.png (both views write the same name; fig_gamma_fine copies it to leiden_resolution_fine_<level>.png)
Notes:  the four-panel view adds cross-seed stability and modularity on log axes. The county compact view fixes the gamma axis to 0.15-0.30 and the cluster axis to 0-2000 to match the sweep grid in _master.do; the tract view derives its limits from the data.
"""

import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter, NullFormatter

SELECTOR_C = "#1f6f3c"
CHOSEN_C = "crimson"

def _at(xs, ys, x0):
    xs = list(xs)
    ys = list(ys)
    i = min(range(len(xs)), key=lambda j: abs(xs[j] - x0))
    return xs[i], ys[i]

def _compact(tables_dir, figures_dir, level, picked_gamma=None):
    df = pd.read_csv(f"{tables_dir}/leiden_resolution_sweep_{level}.csv")
    df = df.sort_values("resolution")
    g = df["resolution"]

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(12, 4.6))

    axL.plot(g, df["n_clusters_est_national"], marker="o", ms=4, lw=1.8,
             color=SELECTOR_C)
    axL.set_xlabel("Leiden resolution γ (CPM)")
    axL.set_ylabel("estimated national #clusters")

    axR.plot(g, df["mean_largest_frac"], marker="o", ms=3, lw=1.6,
             color="firebrick", label="largest-community fraction")
    axR.plot(g, df["mean_singleton_frac"], marker="s", ms=3, lw=1.6,
             color="steelblue", label="singleton fraction")
    axR.axhline(0.50, color="0.5", ls=":", lw=1.2, label="0.50 (degenerate beyond)")
    axR.set_xlabel("Leiden resolution γ (CPM)")
    axR.set_ylabel("fraction of nodes")

    if picked_gamma is not None:
        for ax_ in (axL, axR):
            ax_.axvline(picked_gamma, color=CHOSEN_C, lw=1.5)

        gx, ny = _at(g, df["n_clusters_est_national"], picked_gamma)
        axL.plot([gx], [ny], "o", ms=9, mfc="none", mec=CHOSEN_C, mew=1.6, zorder=5)
        axL.plot([], [], color=CHOSEN_C, lw=1.5, label=f"chosen γ={picked_gamma:g}")
        axL.annotate(f"at chosen: {ny:,.0f}", xy=(gx, ny), xytext=(12, -18),
                     textcoords="offset points", color=CHOSEN_C, fontsize=9,
                     ha="left", va="top")
        axL.legend(fontsize=8, loc="upper left")

        _, blob = _at(g, df["mean_largest_frac"], picked_gamma)
        gx2, sing = _at(g, df["mean_singleton_frac"], picked_gamma)
        axR.plot([gx2, gx2], [blob, sing], "o", ms=9, mfc="none", mec=CHOSEN_C,
                 mew=1.6, zorder=5)
        axR.annotate(f"at chosen γ={picked_gamma:g}:\n"
                     f"singletons {sing:.1%}\nlargest {blob:.1%}",
                     xy=(0.035, 0.64), xycoords="axes fraction", color=CHOSEN_C,
                     fontsize=9, va="top", ha="left")

    axR.legend(fontsize=8, loc="upper right")

    plain = FuncFormatter(lambda v, _: f"{v:g}")
    if level == "county":
        gticks = [0.15, 0.20, 0.25, 0.30]
        for ax_ in (axL, axR):
            ax_.set_xticks(gticks)
            ax_.xaxis.set_major_formatter(plain)
            ax_.xaxis.set_minor_formatter(NullFormatter())
            ax_.set_xlim(0.15, 0.30)
        axL.set_yticks([0, 500, 1000, 1500, 2000])
        axL.set_ylim(top=2000)
        axL.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:,.0f}"))
        axL.yaxis.set_minor_formatter(NullFormatter())
        axR.set_ylim(top=0.7)
    else:
        gmin, gmax = float(g.min()), float(g.max())
        if picked_gamma is not None:
            gmin, gmax = min(gmin, picked_gamma), max(gmax, picked_gamma)
        pad = 0.02 * ((gmax - gmin) or 1.0)
        gticks = np.linspace(gmin, gmax, 6)
        for ax_ in (axL, axR):
            ax_.set_xticks(gticks)
            ax_.xaxis.set_major_formatter(plain)
            ax_.xaxis.set_minor_formatter(NullFormatter())
            ax_.set_xlim(gmin - pad, gmax + pad)

        cmin = max(1.0, float(df["n_clusters_est_national"].min()))
        cmax = float(df["n_clusters_est_national"].max())
        axL.set_yscale("log")
        axL.set_ylim(cmin / 1.5, cmax * 1.5)
        axL.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:,.0f}"))
        axL.yaxis.set_minor_formatter(NullFormatter())

        axR.set_ylim(0, 1.02)

    for ax_ in (axL, axR):
        ax_.grid(alpha=0.25, which="both")

    fig.tight_layout()
    out = f"{figures_dir}/leiden_resolution_lines_{level}.png"
    fig.savefig(out, dpi=140)
    plt.close(fig)
    print("  ->", out)

def main(tables_dir, figures_dir, level, picked_gamma=None):
    df = pd.read_csv(f"{tables_dir}/leiden_resolution_sweep_{level}.csv")
    df = df.sort_values("resolution")
    g = df["resolution"]

    fig, ax = plt.subplots(2, 2, figsize=(13, 9))

    a = ax[0, 0]
    a.loglog(g, df["n_clusters_est_national"], "o-", color="navy")
    a.set_ylabel("estimated national #clusters")
    a.set_title("Granularity vs γ")

    b = ax[0, 1]
    b.semilogx(g, df["mean_stability"], "o-", color="seagreen")
    b.set_ylabel("mean cross-seed stability (NMI)")
    b.set_title("Reproducibility across random seeds")

    c = ax[1, 0]
    c.semilogx(g, df["mean_modularity"], "o-", color="purple")
    c.axhline(0, color="0.6", lw=0.8)
    c.set_ylabel("mean modularity")
    c.set_xlabel("Leiden resolution γ (CPM)")
    c.set_title("Modularity of the partition")

    d = ax[1, 1]
    d.semilogx(g, df["mean_largest_frac"], "o-", color="firebrick",
               label="largest-community fraction")
    d.semilogx(g, df["mean_singleton_frac"], "s-", color="steelblue",
               label="singleton fraction")
    d.axhline(0.50, color="0.5", ls=":", lw=1.2,
              label="0.50 (degenerate beyond)")
    d.set_ylabel("fraction of nodes")
    d.set_xlabel("Leiden resolution γ (CPM)")
    d.set_title("Degeneracy guards")
    d.legend(fontsize=8)

    for ax_ in ax.ravel():
        ax_.grid(alpha=0.25, which="both")
        if picked_gamma is not None:
            ax_.axvline(picked_gamma, color="darkorange", ls="--", lw=1.2)
    if picked_gamma is not None:
        ax[0, 0].annotate(f"picked γ = {picked_gamma:g}",
                          xy=(picked_gamma, 1), xycoords=("data", "axes fraction"),
                          xytext=(4, -12), textcoords="offset points",
                          color="darkorange", fontsize=9)

    fig.suptitle(f"Leiden CPM resolution sweep — {level}", fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    out = f"{figures_dir}/leiden_resolution_lines_{level}.png"
    fig.savefig(out, dpi=140)
    plt.close(fig)
    print("  ->", out)

if __name__ == "__main__":
    args = sys.argv[1:]

    def _opt(flag, default=None, cast=str):
        if flag in args:
            i = args.index(flag)
            v = cast(args[i + 1])
            del args[i:i + 2]
            return v
        return default

    picked = _opt("--picked-gamma", None, float)
    compact = "--compact" in args
    if compact:
        args.remove("--compact")

    if len(args) != 3:
        sys.exit("usage: leiden_resolution_lines.py <tables_dir> <figures_dir> "
                 "<level> [--picked-gamma G] [--compact]")
    render = _compact if compact else main
    render(args[0], args[1], args[2], picked_gamma=picked)
