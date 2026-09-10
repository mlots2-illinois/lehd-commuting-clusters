"""
robinson_figures.py -- figures for the unit-of-analysis (Robinson / MAUP) results of assess_11_robinson.do.

Called from assess/assess_11_robinson.do as:
    python robinson_figures.py --ladder $tables/robinson_ladder_<run_tag>.csv --specs $tables/robinson_specs_<run_tag>.csv --outdir $figures
and again from assess/assess_12_random_null.do with --null $tables/robinson_null_<run_tag>.csv added.
Reads:  --ladder (rel, rellab, unit, family, n_units, b, se, r2, sign, lowpower); --specs (rel, tag, method, n_units, b, se, sign); --null (rel, nullkind, k_target, n_units, b), optional
Writes: in --outdir:
        robinson_heatmap.png and robinson_heatmap_r2.png  relationship x unit bubbles, colour = signed t, size = estimate within row or R2
        robinson_ladder.png                                estimate against number of units for earn_rent, administrative versus functional
        robinson_speccurve.png                             specification curve over every delineation plus the administrative ladder
        robinson_fragility.png                             share of delineations by sign category, per relationship
        robinson_params.png                                estimate against each perturbed calibration parameter, HC and Leiden rows
        robinson_null.png                                  real partitions against the random nulls (needs --null)
        The last four need --specs; without it only the heatmaps and the ladder are written.
Notes:  colour encodes the signed t-statistic or the sign category, never alone: every cell and marker also carries the printed estimate or a glyph, so the figures survive greyscale printing and colour-vision deficiency.
        The sign label "null" is a category, so the CSVs are read with keep_default_na=False. STAR = earn_rent is the relationship carried through in detail.
        The heatmap sizes bubbles within each row because slopes differ by three orders of magnitude across relationships.
"""

import os
import argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
import matplotlib.patheffects as pe

# ColorBrewer PRGn/RdBu endpoints: distinguishable under deuteranopia
POS, NUL, NEG, MISS = "#1B7837", "#B4B4B4", "#B2182B", "#FFFFFF"
COL = {"pos": POS, "null": NUL, "neg": NEG, ".": MISS}
GLYPH = {"pos": "+", "null": "·", "neg": "−", ".": ""}

UNIT_ORDER = ["Census region", "Census division", "State",
              "HC (county)", "Leiden (county)", "County",
              "HC (tract)", "Leiden (tract)", "Tract"]
UNIT_SHORT = ["Region", "Division", "State",
              "HC", "Leiden", "County",
              "HC", "Leiden", "Tract"]
STAR = "earn_rent"          # the relationship carried through in detail


def _relorder(df):
    """Relationships ordered from unit-robust to unit-fragile."""
    piv = df.pivot_table(index="rel", columns="unit", values="sign",
                         aggfunc="first")
    score = {}
    for rel, row in piv.iterrows():
        vals = [v for v in row.dropna().tolist() if v != "."]
        score[rel] = (("pos" in vals and "neg" in vals),
                      len(set(vals)), rel)
    return [r for r in sorted(score, key=lambda r: score[r])]


def _diverging():
    """Red -> pale -> green, using the same endpoints as the categorical fills."""
    from matplotlib.colors import LinearSegmentedColormap
    return LinearSegmentedColormap.from_list(
        "signdiv", [(0.00, NEG), (0.28, "#E8A0A0"), (0.5, "#F0F0F0"),
                    (0.72, "#8CC7A1"), (1.00, POS)])


TCLIP = 12.0          # |t| beyond this is drawn at full saturation


def fig_heatmap(lad, outdir, sizeby="b"):
    """Relationship x unit, as bubbles.

    Two quantities are read at once, so they get separate channels.  COLOUR is
    the signed t-statistic: how confident the estimate is, and in which
    direction.  SIZE is the magnitude of the estimate printed inside the
    bubble, scaled against the largest estimate in the same row, so that a
    reader comparing units within a relationship sees the effect grow and
    shrink directly.  Sizes are not comparable between rows, because the
    relationships are not in the same units.
    """
    lad = lad[lad.unit.isin(UNIT_ORDER)]
    rels = _relorder(lad)
    lab = lad.drop_duplicates("rel").set_index("rel")["rellab"].to_dict()
    # Size encodes the estimate printed inside the bubble.  Slopes differ by
    # three orders of magnitude ACROSS relationships (0.03 to 25,559), so a
    # single absolute scale is impossible; but within a row every estimate is
    # in the same units and directly comparable, which is also the only
    # comparison the figure asks the reader to make.  Each row is therefore
    # scaled by its own largest estimate.
    if sizeby == "r2":
        # "impact on the model" read as share of variance explained.  R-squared
        # is already on a common 0-1 scale, so it needs no per-row rescaling
        # and sizes ARE comparable between rows.
        rowmax = {r: 1.0 for r in lad["rel"].unique()}
    else:
        rowmax = (lad.assign(_a=lad["b"].abs())
                  .groupby("rel")["_a"].max().to_dict())
    n_r, n_u = len(rels), len(UNIT_ORDER)
    cmap = _diverging()
    RMAX = 0.455          # radius of the largest estimate in a row
    RMIN = 0.150          # floor, so the smallest estimate stays legible

    fig, ax = plt.subplots(figsize=(11.8, 0.62 * n_r + 3.0))
    ax.set_facecolor("#FBFBFB")
    for j in range(n_u + 1):
        ax.plot([j, j], [0, n_r], color="#E4E4E4", linewidth=0.7, zorder=1)
    for i in range(n_r + 1):
        ax.plot([0, n_u], [i, i], color="#E4E4E4", linewidth=0.7, zorder=1)

    for i, rel in enumerate(rels):
        denom = rowmax.get(rel)
        for j, unit in enumerate(UNIT_ORDER):
            sub = lad[(lad.rel == rel) & (lad.unit == unit)]
            if sub.empty:
                continue
            r = sub.iloc[0]
            b, se = r["b"], r["se"]
            if not np.isfinite(b) or not np.isfinite(se) or se == 0 or not denom:
                ax.text(j + 0.5, n_r - 1 - i + 0.5, "n/a", ha="center",
                        va="center", fontsize=7, color="#999999", zorder=3)
                continue
            t = b / se
            if sizeby == "r2":
                ratio = float(r["r2"]) if np.isfinite(r.get("r2", np.nan)) else 0.0
            else:
                ratio = abs(b) / denom if denom else 0.0
            rad = RMIN + (RMAX - RMIN) * np.sqrt(min(ratio, 1.0))
            col = cmap(0.5 + 0.5 * np.clip(t / TCLIP, -1, 1))
            low = bool(r.get("lowpower", 0))
            ax.add_patch(plt.Circle((j + 0.5, n_r - 1 - i + 0.5), rad,
                                    facecolor=col,
                                    edgecolor="#404040" if not low else "#404040",
                                    linewidth=0.7,
                                    linestyle="--" if low else "-", zorder=3))
            # the estimate itself, so the figure is still a table
            txt = ("%.0f" % b) if abs(b) >= 10 else ("%.2f" % b)
            lum = 0.299 * col[0] + 0.587 * col[1] + 0.114 * col[2]
            fs = 6.6 if rad >= 0.21 else (5.6 if rad >= 0.175 else 4.9)
            ax.text(j + 0.5, n_r - 1 - i + 0.5, txt, ha="center", va="center",
                    fontsize=fs, zorder=4,
                    color="white" if lum < 0.55 else "#222222")

    ax.set_xlim(0, n_u)
    ax.set_ylim(0, n_r)
    ax.set_xticks(np.arange(n_u) + 0.5)
    ax.set_xticklabels(UNIT_SHORT, fontsize=9)
    ax.set_yticks(np.arange(n_r) + 0.5)
    ax.set_yticklabels([lab.get(r, r) for r in reversed(rels)], fontsize=9)
    ax.tick_params(length=0)
    ax.set_aspect("equal")
    for sp in ax.spines.values():
        sp.set_visible(False)

    # Only the multi-column families get a bracket; the two single-column
    # groups (County, Tract) are already named by their own column label.
    for x0, x1, name, lift in [(0, 3, "Administrative", 0.30),
                               (3, 5, "Functional (county-built)", 0.30),
                               (6, 8, "Functional (tract-built)", 0.30)]:
        ax.plot([x0 + 0.06, x1 - 0.06], [n_r + 0.16] * 2, color="#777777",
                linewidth=1.1, clip_on=False)
        ax.text((x0 + x1) / 2, n_r + lift, name, ha="center", va="bottom",
                fontsize=7.4, color="#444444")

    ax.set_title("Does the areal unit change the conclusion?", fontsize=13,
                 loc="left", pad=62)
    if sizeby == "r2":
        sizenote = ("Size: $R^2$, the share of variance the single regressor "
                    "explains. Sizes are comparable everywhere.")
    else:
        sizenote = ("Size: the estimate itself, scaled within each row. "
                    "Sizes compare across a row, not down a column.")
    ax.text(0, 1.055,
            "Colour: signed t-statistic (confidence and direction). "
            + sizenote
            + "\nState-clustered standard errors; a dashed outline marks "
              "fewer than 30 units.",
            transform=ax.transAxes, fontsize=8.4, color="#555555",
            va="bottom")

    sm = plt.cm.ScalarMappable(cmap=cmap,
                               norm=plt.Normalize(-TCLIP, TCLIP))
    cb = fig.colorbar(sm, ax=ax, orientation="horizontal", fraction=0.033,
                      pad=0.20, aspect=46)
    cb.set_label("t-statistic   (|t| > 1.96 is significant at 5%)",
                 fontsize=8.4)
    cb.ax.tick_params(labelsize=7.6)
    for v in (-1.96, 1.96):
        cb.ax.axvline(v, color="#222222", linewidth=1.0, linestyle=":")

    handles = [Line2D([], [], marker="o", linestyle="none",
                      markerfacecolor="#BBBBBB", markeredgecolor="#404040",
                      markersize=(RMIN + (RMAX - RMIN) * np.sqrt(r)) * 30.0,
                      label=({0.1: "$R^2$ = 0.1", 0.5: "0.5", 1.0: "1.0"}[r]
                             if sizeby == "r2" else
                             {0.1: "10% of the row's largest", 0.5: "50%",
                              1.0: "the row's largest"}[r]))
               for r in (0.1, 0.5, 1.0)]
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.055),
              ncol=3, frameon=False, fontsize=8.2, handletextpad=1.3,
              columnspacing=3.0)
    out = os.path.join(outdir, "robinson_heatmap%s.png"
                       % ("_r2" if sizeby == "r2" else ""))
    fig.savefig(out, dpi=170, bbox_inches="tight")
    plt.close(fig)
    print("  -> %s" % out)


def fig_ladder(lad, spec, outdir):
    """Estimate against number of units, administrative versus functional.

    Two panels because the coarse administrative estimates are an order of
    magnitude away from everything else: the left panel shows the whole range,
    the right zooms on the granularities where administrative and functional
    units can actually be compared like with like.
    """
    a = lad[(lad.rel == STAR) & (lad.family.isin(["Administrative", "Baseline"]))]
    a = a[(a.n_units >= 2) & a.b.notna()].copy()
    fc = lad[(lad.rel == STAR) & (lad.family == "Functional")
             & lad.unit.str.contains(r"\(county\)") & lad.b.notna()]
    f = spec[(spec.rel == STAR) & spec.b.notna()].copy() if spec is not None \
        else pd.DataFrame(columns=["n_units", "b", "se", "method"])
    total = len(a) + len(fc) + len(f)

    fig, (ax, az) = plt.subplots(1, 2, figsize=(13.6, 5.8),
                                 gridspec_kw={"width_ratios": [1.0, 1.15]})

    def draw(target, annotate):
        target.axhline(0, color="#333333", linewidth=1.0, zorder=1)
        if not f.empty:
            target.vlines(f.n_units, f.b - 1.96 * f.se, f.b + 1.96 * f.se,
                          color="#B8B8B8", linewidth=0.7, alpha=0.55, zorder=2)
            for meth, mk, cl, nm in [("hc", "o", "#1F4E79", "HC"),
                                     ("leiden", "^", "#7A2E2E", "Leiden")]:
                g = f[f.method == meth]
                target.scatter(g.n_units, g.b, s=30, marker=mk, alpha=0.8,
                               edgecolor="white", linewidth=0.5, color=cl,
                               zorder=3,
                               label="Functional, %s (%d)" % (nm, len(g)))
        if not fc.empty:
            target.scatter(fc.n_units, fc.b, s=95, marker="D", color="#2E7D32",
                           edgecolor="white", linewidth=1.1, zorder=4,
                           label="Functional, built from counties (%d)" % len(fc))
        target.errorbar(a.n_units, a.b, yerr=1.96 * a.se, fmt="s", ms=10,
                        color="#B2182B", ecolor="#B2182B", elinewidth=1.5,
                        capsize=4, zorder=5, linestyle="none",
                        label="Administrative (%d)" % len(a))
        if annotate:
            for _, r in a.iterrows():
                target.annotate(r["unit"], (r["n_units"], r["b"]),
                                textcoords="offset points", xytext=(0, 15),
                                ha="center", fontsize=8.6, color="#B2182B",
                                weight="bold")
        target.set_xscale("log")
        target.set_xlabel("Number of areal units (log scale)")
        target.grid(axis="y", color="#EDEDED", linewidth=0.8)
        target.set_axisbelow(True)
        for sp in ("top", "right"):
            target.spines[sp].set_visible(False)

    draw(ax, True)
    ax.set_ylabel("Slope: $ mean earnings per +1pp renter-occupied")
    ax.set_title("A.  The whole ladder", fontsize=10.5, loc="left")

    draw(az, False)
    az.set_xlim(400, 1.3e5)
    lo = min(f.b.min() if not f.empty else 0, a.b.min(), fc.b.min()) - 120
    az.set_ylim(lo, 260)
    az.set_title("B.  Where the two families overlap", fontsize=10.5, loc="left")
    az.legend(loc="lower right", frameon=False, fontsize=8.4)

    # the matched-granularity comparison, called out explicitly
    cty = a[a.unit == "County"]
    if not cty.empty and not f.empty:
        cn, cb = float(cty.iloc[0]["n_units"]), float(cty.iloc[0]["b"])
        near = f[(f.n_units > cn * 0.85) & (f.n_units < cn * 1.15)]
        az.annotate("County", (cn, cb), textcoords="offset points",
                    xytext=(0, 16), ha="center", fontsize=9,
                    color="#B2182B", weight="bold")
        if not near.empty:
            nb = float(near.b.mean())
            az.annotate("", xy=(cn, nb), xytext=(cn, cb),
                        arrowprops=dict(arrowstyle="<->", color="#333333",
                                        linewidth=1.3, shrinkA=2, shrinkB=2))
            az.text(cn * 1.1, (cb + nb) / 2,
                    "same number of units,\n%.0f%% of the tract estimate lost"
                    % (100 * (1 - abs(cb) / abs(float(
                        a[a.unit == "Tract"].iloc[0]["b"])))),
                    fontsize=8.2, va="center", color="#333333")

    fig.suptitle("The same relationship, the same records, %d ways of drawing "
                 "the map" % total, fontsize=12.5, x=0.006, ha="left")
    fig.text(0.006, 0.925, "Mean earnings on percent renter-occupied. "
             "Positive above the line, negative below. Bars are 95% "
             "confidence intervals on state-clustered standard errors.",
             fontsize=8.8, color="#555555", ha="left")
    fig.tight_layout(rect=[0, 0, 1, 0.90])
    out = os.path.join(outdir, "robinson_ladder.png")
    fig.savefig(out, dpi=170)
    plt.close(fig)
    print("  -> %s" % out)


def fig_speccurve(spec, lad, outdir):
    f = spec[spec.rel == STAR].copy()
    if f.empty:
        print("  (skip speccurve: no specification rows)")
        return
    extra = lad[(lad.rel == STAR) & (lad.unit.isin(UNIT_ORDER))].copy()
    extra["method"] = "admin"
    extra["tag"] = extra["unit"]
    f = pd.concat([f[["rel", "tag", "method", "n_units", "b", "se", "sign"]],
                   extra[["rel", "tag", "method", "n_units", "b", "se", "sign"]]],
                  ignore_index=True)
    f = f.dropna(subset=["b"]).sort_values("b").reset_index(drop=True)
    x = np.arange(len(f))

    fig, (ax, ax2) = plt.subplots(
        2, 1, figsize=(11.5, 7.4), sharex=True,
        gridspec_kw={"height_ratios": [3.0, 1.0], "hspace": 0.06})
    ax.axhline(0, color="#333333", linewidth=1.0)
    ax.vlines(x, f.b - 1.96 * f.se, f.b + 1.96 * f.se,
              color="#C8C8C8", linewidth=1.0, zorder=2)
    ax.scatter(x, f.b, s=13, c=[COL.get(s, MISS) for s in f["sign"]],
               edgecolor="#555555", linewidth=0.25, zorder=3)
    ax.set_ylabel("Slope (95% CI)")
    ax.set_title("Specification curve: mean earnings on percent "
                 "renter-occupied\nEvery delineation the pipeline produced, "
                 "plus the administrative ladder", fontsize=11.5, loc="left")
    ax.grid(axis="y", color="#EEEEEE", linewidth=0.8)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.legend(handles=[Patch(facecolor=POS, label="positive, significant"),
                       Patch(facecolor=NUL, label="not significant"),
                       Patch(facecolor=NEG, label="negative, significant")],
              loc="upper left", frameon=False, fontsize=8.4)

    rows = [("Administrative unit", f.method == "admin"),
            ("Functional, HC", f.method == "hc"),
            ("Functional, Leiden", f.method == "leiden"),
            ("Coarser than the county", f.n_units < 3100),
            ("Finer than the county", f.n_units >= 3100)]
    for i, (name, mask) in enumerate(rows):
        y = len(rows) - 1 - i
        ax2.scatter(x[mask.values], np.full(mask.sum(), y), s=6,
                    color="#33445A", marker="|")
    ax2.set_yticks(range(len(rows)))
    ax2.set_yticklabels([r[0] for r in reversed(rows)], fontsize=8.2)
    ax2.set_ylim(-0.6, len(rows) - 0.4)
    ax2.set_xlabel("Specifications, ordered by the size of the estimate")
    ax2.tick_params(length=0)
    for s in ("top", "right", "left"):
        ax2.spines[s].set_visible(False)
    fig.tight_layout()
    out = os.path.join(outdir, "robinson_speccurve.png")
    fig.savefig(out, dpi=170)
    plt.close(fig)
    print("  -> %s" % out)


def fig_fragility(spec, lad, outdir):
    if spec is None or spec.empty:
        print("  (skip fragility: no specification rows)")
        return
    rels = _relorder(lad[lad.unit.isin(UNIT_ORDER)])
    lab = lad.drop_duplicates("rel").set_index("rel")["rellab"].to_dict()
    fig, ax = plt.subplots(figsize=(10.2, 0.46 * len(rels) + 2.0))
    for i, rel in enumerate(rels):
        g = spec[spec.rel == rel]
        n = len(g)
        if not n:
            continue
        y = len(rels) - 1 - i
        left = 0.0
        for key in ("neg", "null", "pos"):
            w = 100.0 * (g["sign"] == key).sum() / n
            if w <= 0:
                continue
            ax.barh(y, w, left=left, height=0.62, color=COL[key],
                    edgecolor="white", linewidth=0.8)
            if w >= 8:
                ax.text(left + w / 2, y, "%.0f%%" % w, ha="center",
                        va="center", fontsize=7.6,
                        color="white" if key != "null" else "#444444")
            left += w
    ax.set_yticks(range(len(rels)))
    ax.set_yticklabels([lab.get(r, r) for r in reversed(rels)], fontsize=8.6)
    ax.set_xlim(0, 100)
    ax.set_xlabel("Share of delineations (%)")
    nspec = spec.groupby("rel").size().max()
    ax.set_title("Across %d delineations, the sign almost never changes"
                 % nspec, fontsize=12.5, loc="left", pad=26)
    ax.text(0, 1.008, "Every calibration setting, threshold, temporal window "
            "and division variant the pipeline produced, spanning $\\alpha$ "
            "from 0 (distance only) to 1 (flows only).",
            transform=ax.transAxes, fontsize=8.6, color="#555555", va="bottom")
    ax.tick_params(length=0)
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.legend(handles=[Patch(facecolor=NEG, label="negative, significant"),
                       Patch(facecolor=NUL, label="not significant"),
                       Patch(facecolor=POS, label="positive, significant")],
              loc="upper left", bbox_to_anchor=(0, -0.14), ncol=3,
              frameon=False, fontsize=8.2)
    fig.tight_layout()
    out = os.path.join(outdir, "robinson_fragility.png")
    fig.savefig(out, dpi=170, bbox_inches="tight")
    plt.close(fig)
    print("  -> %s" % out)


# alpha is a weight on [0, 1]; it is swept on its own scale, not as a
# multiple of the calibrated 0.8.  The other three have no intrinsic ceiling
# and are swept as multiples.
ALPHA_TAGS = {"sens_a_full000": 0.0, "sens_a_full020": 0.2,
              "sens_a_m050": 0.4, "sens_a_m075": 0.6,
              "baseline": 0.8, "sens_a_m125": 1.0}
MULT_TAGS = {"m050": 0.50, "m075": 0.75, "m125": 1.25, "m150": 1.50}

# (row method, panel) -> (stem, symbol, plain name).  k enters only the
# hierarchical construction and gamma only Leiden, so each method's row shows
# the two shared parameters plus the one that is its own.
PANELS = {
    "hc":     [("alpha", r"$\alpha$", "flow weight"),
               ("sens_cap", r"$\bar{d}$", "distance cap"),
               ("sens_k", r"$k$", "neighbourhood size")],
    "leiden": [("alpha", r"$\alpha$", "flow weight"),
               ("sens_cap", r"$\bar{d}$", "distance cap"),
               ("sens_res", r"$\gamma$", "resolution")],
}


def fig_params(spec, lad, outdir, rels=None):
    """How the conclusion moves as each calibrated parameter is perturbed."""
    if spec is None or spec.empty:
        print("  (skip params: no specification rows)")
        return
    if rels is None:
        ordered = _relorder(lad[lad.unit.isin(UNIT_ORDER)])
        rels = ordered[-4:][::-1]
    lab = lad.drop_duplicates("rel").set_index("rel")["rellab"].to_dict()
    base_tract = (lad[lad.unit == "Tract"].drop_duplicates("rel")
                  .set_index("rel")["b"].to_dict())
    cmap = plt.cm.viridis
    colors = {r: cmap(i / max(len(rels) - 1, 1)) for i, r in enumerate(rels)}

    fig, axes = plt.subplots(2, 3, figsize=(13.6, 7.4), sharey=True)
    for row, meth in enumerate(("hc", "leiden")):
        for col, (stem, sym, plain) in enumerate(PANELS[meth]):
            ax = axes[row][col]
            g_all = spec[spec.method == meth]
            for rel in rels:
                g = g_all[g_all.rel == rel]
                denom = base_tract.get(rel)
                if not denom:
                    continue
                pts = []
                if stem == "alpha":
                    for tag, val in ALPHA_TAGS.items():
                        rr = g[g.tag == tag]
                        if not rr.empty:
                            pts.append((val, rr.iloc[0]["b"], rr.iloc[0]["se"]))
                else:
                    rr = g[g.tag == "baseline"]
                    if not rr.empty:
                        pts.append((1.0, rr.iloc[0]["b"], rr.iloc[0]["se"]))
                    for suf, m in MULT_TAGS.items():
                        rr = g[g.tag == "%s_%s" % (stem, suf)]
                        if not rr.empty:
                            pts.append((m, rr.iloc[0]["b"], rr.iloc[0]["se"]))
                if len(pts) < 2:
                    continue
                pts.sort()
                xs = [p[0] for p in pts]
                ys = [p[1] / denom for p in pts]
                es = [1.96 * p[2] / abs(denom) for p in pts]
                ax.errorbar(xs, ys, yerr=es, marker="o", ms=4.5, linewidth=1.5,
                            capsize=2.5, color=colors[rel], alpha=0.92,
                            label=lab.get(rel, rel) if (row == 0 and col == 0)
                            else None)
            ax.axhspan(-0.07, 0.07, color="#F6DADA", zorder=0)
            ax.axhline(0, color="#B2182B", linewidth=1.1, zorder=1)
            ax.axhline(1, color="#333333", linewidth=1.0, linestyle="--",
                       zorder=1)
            if stem == "alpha":
                ax.set_xlim(-0.05, 1.05)
                ax.set_xticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
                ax.axvline(0.8, color="#AAAAAA", linewidth=0.9, linestyle=":")
                ax.set_xlabel("$\\alpha$ (calibrated: 0.8)", fontsize=9)
            else:
                ax.set_xlim(0.42, 1.58)
                ax.set_xticks([0.5, 0.75, 1.0, 1.25, 1.5])
                ax.axvline(1.0, color="#AAAAAA", linewidth=0.9, linestyle=":")
                ax.set_xlabel("multiple of the calibrated value", fontsize=9)
            ax.set_title("%s  %s" % (sym, plain), fontsize=10.5)
            ax.grid(axis="y", color="#EFEFEF", linewidth=0.8)
            ax.set_axisbelow(True)
            for sp in ("top", "right"):
                ax.spines[sp].set_visible(False)
        axes[row][0].set_ylabel(("Hierarchical clustering" if meth == "hc"
                                 else "Leiden")
                                + "\nestimate / tract estimate", fontsize=10)
    axes[0][0].text(0.02, 1.0, " reproduces the tract answer", fontsize=7.4,
                    color="#333333", va="bottom")
    axes[0][0].text(0.02, 0.0, " sign flip", fontsize=7.4, color="#B2182B",
                    va="bottom")
    h, l = axes[0][0].get_legend_handles_labels()
    fig.legend(handles=h, labels=l, loc="lower center", ncol=min(len(l), 4),
               frameon=False, fontsize=8.6, bbox_to_anchor=(0.5, -0.055))
    fig.suptitle("Does the conclusion survive the calibration?", fontsize=13,
                 x=0.006, ha="left")
    fig.text(0.006, 0.945, "Each calibrated parameter perturbed one at a time, "
             "for the four relationships most sensitive to the areal unit. "
             "Rows separate the two algorithms; 1.0 reproduces the "
             "tract-level answer and 0 is a sign flip.",
             fontsize=8.8, color="#555555", ha="left")
    fig.tight_layout(rect=[0, 0, 1, 0.925])
    out = os.path.join(outdir, "robinson_params.png")
    fig.savefig(out, dpi=170, bbox_inches="tight")
    plt.close(fig)
    print("  -> %s" % out)


def fig_null(lad, spec, nul, outdir):
    """Where the real partitions sit against partitions carrying no information.

    The contiguous null is the one that matters: it is spatially connected and
    compact like a real geography, but its boundaries are drawn at random
    rather than from commuting.  A delineation that merely subdivided space
    would land inside that band.
    """
    if nul is None or nul.empty:
        print("  (skip null: no null rows)")
        return
    rel = STAR
    denom = float(lad[(lad.rel == rel) & (lad.unit == "Tract")].b.iloc[0])
    lab = lad[lad.rel == rel].rellab.iloc[0]

    fig, ax = plt.subplots(figsize=(11.2, 5.9))
    ax.axhline(1.0, color="#333333", linestyle="--", linewidth=1.1, zorder=2)
    ax.axhline(0.0, color="#B2182B", linewidth=1.1, zorder=2)

    styles = {"unconstrained": ("#8C8C8C", "Random, unconstrained"),
              "contiguous":    ("#E08214", "Random, spatially contiguous")}
    for kind, (col, name) in styles.items():
        g = nul[(nul.rel == rel) & (nul.nullkind == kind)]
        if g.empty:
            continue
        agg = (g.groupby("k_target")
                 .agg(n=("n_units", "mean"), m=("b", "mean"), sd=("b", "std"))
                 .reset_index().sort_values("n"))
        lo = (agg.m - 2 * agg.sd) / denom
        hi = (agg.m + 2 * agg.sd) / denom
        ax.fill_between(agg.n, lo, hi, color=col, alpha=0.22, zorder=1,
                        label="%s (mean $\\pm$ 2 sd, %d draws)"
                              % (name, int(g.groupby("k_target").size().mean())))
        ax.plot(agg.n, agg.m / denom, color=col, linewidth=1.8, marker="s",
                ms=6, zorder=3)

    f = spec[(spec.rel == rel) & spec.b.notna()]
    ax.scatter(f.n_units, f.b / denom, s=30, marker="o", color="#1F4E79",
               alpha=0.8, edgecolor="white", linewidth=0.5, zorder=4,
               label="Functional delineations (%d)" % len(f))
    fc = lad[(lad.rel == rel) & (lad.family == "Functional")
             & lad.unit.str.contains(r"\(county\)") & lad.b.notna()]
    if not fc.empty:
        ax.scatter(fc.n_units, fc.b / denom, s=90, marker="D",
                   color="#2E7D32", edgecolor="white", linewidth=1.0, zorder=5,
                   label="Functional, built from counties")
    a = lad[(lad.rel == rel) & (lad.family.isin(["Administrative", "Baseline"]))
            & lad.b.notna() & (lad.n_units >= 2)]
    ax.scatter(a.n_units, a.b / denom, s=110, marker="s", color="#B2182B",
               edgecolor="white", linewidth=1.0, zorder=6,
               label="Administrative")
    for _, r in a.iterrows():
        ax.annotate(r["unit"], (r["n_units"], r["b"] / denom),
                    textcoords="offset points", xytext=(0, 14), ha="center",
                    fontsize=8.6, color="#B2182B", weight="bold")

    ax.set_xscale("log")
    ax.set_xlabel("Number of areal units (log scale)")
    ax.set_ylabel("Estimate as a fraction of the tract-level estimate")
    ax.set_ylim(-0.35, 1.45)
    ax.set_title("Would any partition have done as well?", fontsize=13,
                 loc="left", pad=30)
    ax.text(0, 1.015, "%s. 1.0 reproduces the tract answer; 0 is a sign flip. "
            "Region, division and state fall far below the axis (their "
            "estimates carry the opposite sign) and are omitted here; see the "
            "ladder figure." % lab,
            transform=ax.transAxes, fontsize=8.6, color="#555555",
            va="bottom", wrap=True)
    ax.grid(axis="y", color="#EFEFEF", linewidth=0.8)
    ax.set_axisbelow(True)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, -0.30),
              ncol=2, frameon=False, fontsize=8.4)
    fig.tight_layout()
    out = os.path.join(outdir, "robinson_null.png")
    fig.savefig(out, dpi=170, bbox_inches="tight")
    plt.close(fig)
    print("  -> %s" % out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ladder", required=True)
    ap.add_argument("--specs", default="")
    ap.add_argument("--null", default="")
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    # keep_default_na=False matters: the sign category "null" is a LABEL, and
    # pandas would otherwise read it as a missing value and blank the cell.
    lad = pd.read_csv(args.ladder, keep_default_na=False, na_values=[""])
    lad["sign"] = lad["sign"].fillna(".").astype(str)
    if "lowpower" not in lad.columns:
        lad["lowpower"] = 0
    lad["lowpower"] = lad["lowpower"].fillna(0)
    spec = None
    if args.specs and os.path.exists(args.specs):
        spec = pd.read_csv(args.specs, keep_default_na=False,
                           na_values=[""])
        spec["sign"] = spec["sign"].fillna(".").astype(str)
    nul = None
    if args.null and os.path.exists(args.null):
        nul = pd.read_csv(args.null, keep_default_na=False, na_values=[""])
    if not os.path.isdir(args.outdir):
        os.makedirs(args.outdir)

    fig_heatmap(lad, args.outdir)
    fig_heatmap(lad, args.outdir, sizeby="r2")
    fig_ladder(lad, spec, args.outdir)
    if spec is not None:
        fig_speccurve(spec, lad, args.outdir)
        fig_fragility(spec, lad, args.outdir)
        fig_params(spec, lad, args.outdir)
        fig_null(lad, spec, nul, args.outdir)


if __name__ == "__main__":
    main()
