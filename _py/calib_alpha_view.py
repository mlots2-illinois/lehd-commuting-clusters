"""
calib_alpha_view.py -- two-panel figure of the alpha choice: self-containment against alpha (selector) and the D dispersion diagnostics at the chosen cap (guard).

Called from calib/calib_01_alpha_cap.do as:
    python calib_alpha_view.py <tables_dir> <figures_dir> <level> <chosen_cap> <chosen_alpha> [<self_containment_csv>]
Reads:  <tables_dir>/calib_dispersion_<level>.csv; the optional self-containment CSV ($tables/partition_quality_alpha_calib_<level>.csv, else partition_quality_alpha_<run_tag>.csv), filtered to geo_level == <level>
Writes: <figures_dir>/calib_alpha_view_<level>.png
Notes:  the dispersion rows are taken at the grid cap nearest the chosen cap. Without the self-containment CSV the left panel carries a placeholder text.
"""

import os
import sys
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SELECTOR_C = "#1f6f3c"
CHOSEN_C = "crimson"

def _mark_at(ax, xs, ys, x0, fmt="{:.3f}"):
    xs = list(xs)
    ys = list(ys)
    if not xs:
        return
    i = min(range(len(xs)), key=lambda j: abs(xs[j] - x0))
    ax.plot([xs[i]], [ys[i]], "o", ms=8, mfc="none", mec=CHOSEN_C, mew=1.6,
            zorder=5)
    ax.annotate(f"at chosen: {fmt.format(ys[i])}", xy=(xs[i], ys[i]),
                xytext=(6, 8), textcoords="offset points", color=CHOSEN_C,
                fontsize=9)

def _dispersion_at_cap(tables_dir, level, chosen_cap):
    df = pd.read_csv(f"{tables_dir}/calib_dispersion_{level}.csv")
    caps = sorted(df["cap_km"].unique())
    cap = min(caps, key=lambda c: abs(c - chosen_cap))
    return df[df["cap_km"] == cap].sort_values("alpha"), cap

def _self_containment(csv_path, level):
    if not csv_path or not os.path.exists(csv_path):
        return None
    sc = pd.read_csv(csv_path)
    sc = sc[sc["geo_level"] == level].sort_values("alpha")
    return sc if len(sc) else None

def plot_alpha_view(tables_dir, figures_dir, level, chosen_cap, chosen_alpha,
                    sc_csv=None):
    disp, cap_used = _dispersion_at_cap(tables_dir, level, chosen_cap)
    if disp.empty:
        print(f"  (skip {level}: no dispersion rows near cap={chosen_cap})")
        return
    sc = _self_containment(sc_csv, level)

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(12, 4.6))

    if sc is not None:
        axL.plot(sc["alpha"], sc["self_containment"], marker="o", ms=4,
                 lw=1.8, color=SELECTOR_C)
        _mark_at(axL, sc["alpha"], sc["self_containment"], chosen_alpha)
        axL.set_ylabel("self-containment (within-cluster flow share)")
    else:
        axL.text(0.5, 0.5,
                 "self-containment not available yet\n"
                 "(run calib_01 Phase D / assess_08)",
                 ha="center", va="center", fontsize=10, color="gray",
                 transform=axL.transAxes)
    axL.axvline(chosen_alpha, color=CHOSEN_C, lw=1.5,
                label=f"chosen α={chosen_alpha:g}")
    axL.set_xlabel("flow weight α")
    axL.legend(fontsize=8, loc="best")
    axL.grid(alpha=0.25)

    axR.plot(disp["alpha"], disp["d_mean"],    marker="o", ms=3, lw=1.6,
             label="mean $D$")
    axR.plot(disp["alpha"], disp["d_sd"],      marker="d", ms=3, lw=1.6,
             label="SD $D$")
    axR.plot(disp["alpha"], disp["d_frac_hi"], marker="s", ms=3, lw=1.6,
             label="share $D>0.95$")
    axR.plot(disp["alpha"], disp["d_frac_lo"], marker="^", ms=3, lw=1.6,
             label="share $D<0.05$")
    axR.axvline(chosen_alpha, color=CHOSEN_C, lw=1.5)
    axR.set_xlabel("flow weight α")
    axR.set_ylabel("$D$ scale diagnostics")
    axR.legend(fontsize=8, loc="best")
    axR.grid(alpha=0.25)

    fig.tight_layout()
    out = f"{figures_dir}/calib_alpha_view_{level}.png"
    fig.savefig(out, dpi=140)
    plt.close(fig)
    print("  ->", out)

if __name__ == "__main__":
    args = sys.argv[1:]
    if len(args) not in (5, 6):
        sys.exit("usage: calib_alpha_view.py <tables_dir> <figures_dir> "
                 "<level> <chosen_cap> <chosen_alpha> [self_containment_csv]")
    tables_dir, figures_dir, level = args[0], args[1], args[2]
    chosen_cap = float(args[3])
    chosen_alpha = float(args[4])
    sc_csv = args[5] if len(args) == 6 else None
    plot_alpha_view(tables_dir, figures_dir, level, chosen_cap, chosen_alpha,
                    sc_csv)
