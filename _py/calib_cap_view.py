"""
calib_cap_view.py -- two-panel figure of the distance-cap choice: P-weighted flow coverage against the cap (selector) and D dispersion at the chosen alpha (guard).

Called from calib/calib_01_alpha_cap.do as:
    python calib_cap_view.py <tables_dir> <figures_dir> <level> <chosen_cap> <chosen_alpha>
Reads:  <tables_dir>/calib_flow_coverage_<level>.csv (from calib_flow_coverage.py); <tables_dir>/calib_dispersion_<level>.csv
Writes: <figures_dir>/calib_cap_view_<level>.png
Notes:  skips the level when the coverage CSV is absent. Dispersion rows are taken at the grid alpha nearest the chosen alpha.
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
                xytext=(6, -14), textcoords="offset points", color=CHOSEN_C,
                fontsize=9)

def _dispersion_at_alpha(tables_dir, level, chosen_alpha):
    df = pd.read_csv(f"{tables_dir}/calib_dispersion_{level}.csv")
    alphas = sorted(df["alpha"].unique())
    a = min(alphas, key=lambda x: abs(x - chosen_alpha))
    return df[df["alpha"] == a].sort_values("cap_km"), a

def plot_cap_view(tables_dir, figures_dir, level, chosen_cap, chosen_alpha):
    cov_csv = f"{tables_dir}/calib_flow_coverage_{level}.csv"
    if not os.path.exists(cov_csv):
        print(f"  (skip {level}: {cov_csv} not found)")
        return
    cov = pd.read_csv(cov_csv).sort_values("cap_km")
    disp, alpha_used = _dispersion_at_alpha(tables_dir, level, chosen_alpha)

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(12, 4.6))

    axL.plot(cov["cap_km"], cov["coverage"], marker="o", ms=4, lw=1.8,
             color=SELECTOR_C, label="pooled coverage")
    _mark_at(axL, cov["cap_km"], cov["coverage"], chosen_cap)
    axL.set_xlabel("distance cap  $\\bar d$  (km)")
    axL.set_ylabel("flow coverage (P-weighted share with $d \\leq \\bar d$)")
    axL.axvline(chosen_cap, color=CHOSEN_C, lw=1.5, label=f"chosen {chosen_cap:g} km")
    axL.legend(fontsize=8, loc="best")
    axL.grid(alpha=0.25)

    axR.plot(disp["cap_km"], disp["d_frac_hi"], marker="s", ms=3, lw=1.6,
             label="share $D>0.95$ (want low)")
    axR.plot(disp["cap_km"], disp["d_sd"], marker="d", ms=3, lw=1.6,
             label="SD $D$ (want large)")
    axR.axvline(chosen_cap, color=CHOSEN_C, lw=1.5)
    axR.set_xlabel("distance cap  $\\bar d$  (km)")
    axR.set_ylabel("$D$ scale diagnostics")
    axR.legend(fontsize=8, loc="best")
    axR.grid(alpha=0.25)

    fig.tight_layout()
    out = f"{figures_dir}/calib_cap_view_{level}.png"
    fig.savefig(out, dpi=140)
    plt.close(fig)
    print("  ->", out)

if __name__ == "__main__":
    args = sys.argv[1:]
    if len(args) != 5:
        sys.exit("usage: calib_cap_view.py <tables_dir> <figures_dir> <level> "
                 "<chosen_cap> <chosen_alpha>")
    plot_cap_view(args[0], args[1], args[2], float(args[3]), float(args[4]))
