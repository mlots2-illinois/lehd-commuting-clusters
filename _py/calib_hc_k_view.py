"""
calib_hc_k_view.py -- four-panel figure of the HC k choice: self-containment (selector), silhouette and modularity (validity), bootstrap stability, and within-cluster spread (guard).

Called from calib/calib_03_hc_k.do and again from calib/calib_05_validity.do as:
    python calib_hc_k_view.py <tables_dir> <figures_dir> <level> <chosen_k>
Reads:  <tables_dir>/partition_quality_k_<level>.csv (required); <tables_dir>/calib_hc_k_summary_<level>.csv and <tables_dir>/hc_stability_<level>.csv (optional)
Writes: <figures_dir>/calib_hc_k_view_<level>.png
Notes:  calib_05 re-runs it after hc_stability.py so the stability panel fills in; the calib_03 call draws that panel as a placeholder.
        Reference lines: silhouette 0.25 (Rousseeuw's no-substantial-structure line) and Jaccard 0.75 (Hennig's stability threshold).
"""

import os
import sys
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SELECTOR_C = "#1f6f3c"
CHOSEN_C = "crimson"
VALID_C = "#2b5d8a"
STAB_C = "#8a5a2b"


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


def plot_hc_k_view(tables_dir, figures_dir, level, chosen_k):
    sc_csv = f"{tables_dir}/partition_quality_k_{level}.csv"
    sp_csv = f"{tables_dir}/calib_hc_k_summary_{level}.csv"
    st_csv = f"{tables_dir}/hc_stability_{level}.csv"
    if not os.path.exists(sc_csv):
        print(f"  (skip {level}: {sc_csv} not found)")
        return
    sc = pd.read_csv(sc_csv)
    sc = sc[sc["geo_level"] == level].sort_values("k_national")
    if sc.empty:
        print(f"  (skip {level}: no self-containment rows)")
        return

    fig, axes = plt.subplots(1, 4, figsize=(20, 4.6))
    axL, axV, axS, axR = axes

    # (1) selector: self-containment
    axL.plot(sc["k_national"], sc["self_containment"], marker="o", ms=4, lw=1.8,
             color=SELECTOR_C)
    _mark_at(axL, sc["k_national"], sc["self_containment"], chosen_k)
    axL.axvline(chosen_k, color=CHOSEN_C, lw=1.5, label=f"chosen k={chosen_k:g}")
    axL.set_xlabel("national k target")
    axL.set_ylabel("self-containment (within-cluster flow share)")
    axL.set_title("selector: self-containment", fontsize=10)
    axL.legend(fontsize=8, loc="best")
    axL.grid(alpha=0.25)

    # (2) internal validity: silhouette (Rousseeuw 1987) + modularity
    has_sil = "silhouette" in sc and sc["silhouette"].notna().any()
    if has_sil:
        axV.plot(sc["k_national"], sc["silhouette"], marker="o", ms=4, lw=1.8,
                 color=VALID_C, label="mean silhouette width")
        _mark_at(axV, sc["k_national"], sc["silhouette"], chosen_k)
        if "modularity" in sc and sc["modularity"].notna().any():
            ax2 = axV.twinx()
            ax2.plot(sc["k_national"], sc["modularity"], marker="s", ms=3,
                     lw=1.4, ls="--", color="gray", label="modularity")
            ax2.set_ylabel("modularity", color="gray", fontsize=9)
            ax2.tick_params(axis="y", labelcolor="gray")
        axV.axhline(0.25, color="gray", lw=1.0, ls=":")
        axV.annotate("0.25: Rousseeuw's 'no substantial structure' line",
                     xy=(0.02, 0.25), xycoords=("axes fraction", "data"),
                     fontsize=7, color="gray", va="bottom")
        axV.axvline(chosen_k, color=CHOSEN_C, lw=1.5)
        axV.set_ylabel("mean silhouette width")
        axV.legend(fontsize=8, loc="best")
    else:
        axV.text(0.5, 0.5, "silhouette not computed\n(run calib_03 with PQ_FULL=1)",
                 ha="center", va="center", fontsize=10, color="gray",
                 transform=axV.transAxes)
    axV.set_xlabel("national k target")
    axV.set_title("internal validity", fontsize=10)
    axV.grid(alpha=0.25)

    # (3) clusterwise bootstrap stability (Hennig 2007)
    if os.path.exists(st_csv):
        st = pd.read_csv(st_csv).sort_values("k_national")
        axS.plot(st["k_national"], st["jaccard_mean"], marker="o", ms=4, lw=1.8,
                 color=STAB_C, label="mean clusterwise Jaccard")
        axS.plot(st["k_national"], st["frac_stable"], marker="^", ms=3, lw=1.3,
                 ls="--", color=STAB_C, alpha=0.6,
                 label="share of clusters $\\geq 0.75$")
        _mark_at(axS, st["k_national"], st["jaccard_mean"], chosen_k)
        axS.axhline(0.75, color="gray", lw=1.0, ls=":")
        axS.annotate("0.75: Hennig's stability threshold",
                     xy=(0.02, 0.75), xycoords=("axes fraction", "data"),
                     fontsize=7, color="gray", va="bottom")
        axS.axvline(chosen_k, color=CHOSEN_C, lw=1.5)
        axS.set_ylim(0, 1)
        axS.legend(fontsize=8, loc="lower right")
    else:
        axS.text(0.5, 0.5, "stability not found\n(run calib_05_validity)",
                 ha="center", va="center", fontsize=10, color="gray",
                 transform=axS.transAxes)
    axS.set_xlabel("national k target")
    axS.set_ylabel("bootstrap stability")
    axS.set_title("clusterwise stability", fontsize=10)
    axS.grid(alpha=0.25)

    # (4) spread guard
    if os.path.exists(sp_csv):
        sp = pd.read_csv(sp_csv).sort_values("k_target_national")
        axR.plot(sp["k_target_national"], sp["mean_bbox"], marker="o", ms=3,
                 lw=1.6, label="mean bbox diagonal (km)")
        axR.plot(sp["k_target_national"], sp["max_bbox"], marker="s", ms=3,
                 lw=1.6, label="max bbox diagonal (km)")
        axR.axvline(chosen_k, color=CHOSEN_C, lw=1.5)
        axR.set_xlabel("national k target")
        axR.set_ylabel("within-cluster spread (km)")
        axR.legend(fontsize=8, loc="best")
    else:
        axR.text(0.5, 0.5, "spread summary not found\n"
                 "(calib_hc_k_summary CSV missing)",
                 ha="center", va="center", fontsize=10, color="gray",
                 transform=axR.transAxes)
    axR.set_title("spread guard", fontsize=10)
    axR.grid(alpha=0.25)

    fig.tight_layout()
    out = f"{figures_dir}/calib_hc_k_view_{level}.png"
    fig.savefig(out, dpi=140)
    plt.close(fig)
    print("  ->", out)


if __name__ == "__main__":
    args = sys.argv[1:]
    if len(args) != 4:
        sys.exit("usage: calib_hc_k_view.py <tables_dir> <figures_dir> <level> "
                 "<chosen_k>")
    plot_hc_k_view(args[0], args[1], args[2], float(args[3]))
