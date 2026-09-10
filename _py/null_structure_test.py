"""
null_structure_test.py -- tests whether the commuting network carries structure beyond geography: observed self-containment against a geography-only cut and a within-division permutation null.

Called from calib/calib_05_validity.do as:
    NT_KGRID=<k1,...> NT_SAMPLE_DIVS=<d1,...> NT_NPERM=<R> python null_structure_test.py <level> <divlookup_dta> <pflow_dta> <cap_km> <alpha> <seed> <out_csv>
    env: NT_KGRID (required), NT_SAMPLE_DIVS (default all), NT_NPERM (default 20; calib_05 passes $nt_nperm = 99), NT_NPROC, PQ_MINFLOW
Reads:  <divlookup_dta>; <pflow_dta> (through partition_quality_alpha)
Writes: <out_csv>, one row per k_national with sc_<prop|raw>_observed, _geo_only, _null_mean, _null_sd, z_<prop|raw>_vs_null, p_<prop|raw>_emp, lift_<prop|raw>_over_geo, n_perm
Notes:  all three partitions are scored on the real flow matrix, so the numbers are comparable to the self-containment calib_03 and assess_08 report.
        geography-only clusters on D with alpha = 0. The permutation null relabels the units of each division's flow matrix at random, so the network keeps its topology and weight distribution and loses only its alignment with geography; repeated R times.
        Numerators and denominators are summed across divisions before the ratio, so the national statistic is a flow share, not a mean of shares. The empirical p is (count(null >= observed) + 1) / (R + 1).
        This is the homogeneity-against-clustering check of Hennig and Lin (2015).
"""

import os
import sys
from collections import defaultdict
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage

sys.path.insert(0, str(Path(__file__).resolve().parent))
import partition_quality_alpha as pq

_W = {}


def _init(div_units, flow_by_div, cap_km, alpha, k_list, total_n, nperm, seed):
    _W.update(div_units=div_units, flow_by_div=flow_by_div, cap_km=cap_km,
              alpha=alpha, k_list=k_list, total_n=total_n, nperm=nperm,
              seed=seed)


def _cut_and_score(d_arr, k_list, n, total_n, i_idx, j_idx, raw_flow, raw_prop):
    """Cluster on d_arr; score every k on the REAL flow matrix.

    Returns (prop_num, prop_den, raw_num, raw_den) per k: the proportional
    score (P_ij, what the dissimilarity is built on) and the raw-volume
    score (f_ij) side by side.
    """
    Z = linkage(d_arr, method="average")
    out = {}
    for k_nat in k_list:
        k_div = max(2, min(int(round(n / (total_n / k_nat))), n))
        lab = fcluster(Z, t=k_div, criterion="maxclust")
        same = lab[i_idx] == lab[j_idx]
        out[k_nat] = (float(raw_prop[same].sum()), float(raw_prop.sum()),
                      float(raw_flow[same].sum()), float(raw_flow.sum()))
    return out


def _process_div(d):
    n_units = len(_W["div_units"][d])
    if n_units < 10:
        return d, n_units, {}

    sub = _W["div_units"][d]
    flows = _W["flow_by_div"].get(d)
    cap_km, alpha = _W["cap_km"], _W["alpha"]
    k_list, total_n = _W["k_list"], _W["total_n"]

    n, i_idx, j_idx, fterm, gterm, raw_flow, raw_prop = \
        pq._division_base(sub, flows, cap_km, d)

    def _score(dd):
        return _cut_and_score(dd, k_list, n, total_n, i_idx, j_idx,
                              raw_flow, raw_prop)

    res = {"observed": _score(alpha * fterm + (1 - alpha) * gterm),
           "geo_only": _score(gterm.copy()),
           "null": []}

    if flows is not None:
        codes = sub["code"].to_numpy()
        ci, cj, P, f = flows
        li = np.searchsorted(codes, ci)
        lj = np.searchsorted(codes, cj)
        rng = np.random.default_rng(_W["seed"] + 7919 * int(d))
        for _ in range(_W["nperm"]):
            perm = rng.permutation(n)
            a = perm[li]
            b = perm[lj]
            lo = np.minimum(a, b)
            hi = np.maximum(a, b)
            rank = pq._triu_rank(lo, hi, n)
            ft = np.ones(len(i_idx), dtype=np.float64)
            ft[rank] = 1.0 - np.minimum(P, 1.0)
            res["null"].append(_score(alpha * ft + (1 - alpha) * gterm))
    return d, n, res


def main(level, divlookup, pflow, cap_km, alpha, seed, out_csv):
    kgrid = os.environ.get("NT_KGRID", "")
    if not kgrid:
        sys.exit("null_structure_test: NT_KGRID not set")
    k_list = [float(x) for x in kgrid.split(",")]
    nperm = int(os.environ.get("NT_NPERM", "20"))
    samp = os.environ.get("NT_SAMPLE_DIVS", "")
    sample_set = {int(x) for x in samp.replace(",", " ").split()} if samp else None

    div_units, code_of, dvA, dvB, total_n = pq._load_divisions(divlookup)
    print(f"  [{level}] {len(div_units)} divisions, {total_n:,} units; "
          f"{len(k_list)} k values; R={nperm} alpha={alpha} cap={cap_km}",
          flush=True)
    flow_by_div = pq._load_flows_by_div(pflow, code_of, dvA, dvB)

    work = [d for d in sorted(div_units)
            if sample_set is None or d in sample_set]
    nproc = max(1, min(int(os.environ.get("NT_NPROC", os.cpu_count() or 1)),
                       len(work)))
    print(f"  [{level}] permuting {len(work)} divisions across {nproc} workers ...",
          flush=True)

    initargs = (div_units, flow_by_div, cap_km, alpha, k_list, total_n,
                nperm, seed)
    if nproc == 1:
        _init(*initargs)
        results = (_process_div(d) for d in work)
        pool = None
    else:
        pool = Pool(processes=nproc, initializer=_init, initargs=initargs)
        results = pool.imap_unordered(_process_div, work)

    # accumulate flow numerators/denominators across divisions, so the
    # national statistic is a genuine flow share rather than a mean of shares
    obs = defaultdict(lambda: np.zeros(4))
    geo = defaultdict(lambda: np.zeros(4))
    nul = defaultdict(lambda: defaultdict(lambda: np.zeros(4)))
    for d, n, res in results:
        if not res:
            continue
        print(f"  div {d}: n={n}", flush=True)
        for k, v in res["observed"].items():
            obs[k] += np.asarray(v)
        for k, v in res["geo_only"].items():
            geo[k] += np.asarray(v)
        for r, rep in enumerate(res["null"]):
            for k, v in rep.items():
                nul[k][r] += np.asarray(v)
    if pool is not None:
        pool.close(); pool.join()

    out = []
    for k in sorted(obs):
        row = {"geo_level": level, "alpha": alpha, "k_national": k}
        # a = 0 -> proportional score (P_ij); a = 2 -> raw volume score (f_ij)
        for tag, a in (("prop", 0), ("raw", 2)):
            o = obs[k][a] / obs[k][a + 1]
            g = geo[k][a] / geo[k][a + 1]
            draws = np.array([v[a] / v[a + 1] for v in nul[k].values()],
                             dtype=float)
            mu, sd = float(draws.mean()), float(draws.std(ddof=1))
            row.update({
                f"sc_{tag}_observed": o, f"sc_{tag}_geo_only": g,
                f"sc_{tag}_null_mean": mu, f"sc_{tag}_null_sd": sd,
                f"z_{tag}_vs_null": (o - mu) / sd if sd > 0 else float("nan"),
                f"p_{tag}_emp": float((draws >= o).sum() + 1) / (len(draws) + 1),
                f"lift_{tag}_over_geo": o - g,
            })
        row["n_perm"] = len(nul[k])
        out.append(row)
    pd.DataFrame(out).to_csv(out_csv, index=False)
    print(f"  wrote {len(out)} rows to {out_csv}", flush=True)


if __name__ == "__main__":
    if len(sys.argv) != 8:
        sys.exit(__doc__)
    main(level=sys.argv[1], divlookup=sys.argv[2], pflow=sys.argv[3],
         cap_km=float(sys.argv[4]), alpha=float(sys.argv[5]),
         seed=int(sys.argv[6]), out_csv=sys.argv[7])
