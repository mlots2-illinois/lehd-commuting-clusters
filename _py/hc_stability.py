"""
hc_stability.py -- clusterwise bootstrap stability of the hierarchical cut (Hennig 2007) across a national k grid.

Called from calib/calib_05_validity.do as:
    HS_KGRID=<k1,...> HS_SAMPLE_DIVS=<d1,...> HS_NBOOT=<B> python hc_stability.py <level> <divlookup_dta> <pflow_dta> <cap_km> <alpha> <seed> <out_csv>
    env: HS_KGRID (required), HS_SAMPLE_DIVS (default all), HS_NBOOT (default 50; calib_05 passes $hs_nboot), HS_SUBFRAC (default 0.75), HS_NPROC, PQ_MINFLOW
Reads:  <divlookup_dta>; <pflow_dta> (through partition_quality_alpha._load_flows_by_div)
Writes: <out_csv>, one row per k_national: jaccard_mean, jaccard_wmean, frac_stable, frac_dissolved, n_clusters, n_div_used, n_units, n_boot, subfrac
Notes:  for each k the division tree is cut at the proportionally allocated k_div that cluster_01 uses, then re-cut on B subsamples of the division's units.
        Each reference cluster is scored by its maximum Jaccard similarity to any cluster of the subsample partition, averaged over replicates (the clusterboot statistic).
        Hennig's conventions: mean Jaccard >= 0.75 is a stable cluster, <= 0.5 a dissolved one.
        This is the hierarchical-arm analogue of the cross-seed reproducibility calib_02 reports for Leiden; it is not monotone in k, so unlike self-containment it can select an interior target.
        The RNG is seeded with <seed> + div_id, so each division is reproducible on its own. Divisions with fewer than 10 units are skipped.
"""

import os
import sys
from collections import defaultdict
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.sparse import coo_matrix
from scipy.spatial.distance import squareform

sys.path.insert(0, str(Path(__file__).resolve().parent))
import partition_quality_alpha as pq


def _max_jaccard(ref, rep):
    """Per-reference-cluster max Jaccard against any replicate cluster."""
    ru, ri = np.unique(ref, return_inverse=True)
    pu, pi = np.unique(rep, return_inverse=True)
    cont = coo_matrix((np.ones(len(ri)), (ri, pi)),
                      shape=(len(ru), len(pu))).toarray()
    a = cont.sum(axis=1, keepdims=True)      # |A_i|
    b = cont.sum(axis=0, keepdims=True)      # |B_j|
    jac = cont / (a + b - cont)
    return ru, jac.max(axis=1)


_W = {}


def _init(div_units, flow_by_div, cap_km, alpha, k_list, total_n,
          nboot, subfrac, seed):
    _W.update(div_units=div_units, flow_by_div=flow_by_div, cap_km=cap_km,
              alpha=alpha, k_list=k_list, total_n=total_n, nboot=nboot,
              subfrac=subfrac, seed=seed)


def _process_div(d):
    n_units = len(_W["div_units"][d])
    if n_units < 10:
        return d, n_units, []

    n, iu, ju, fterm, gterm, _raw, _prop = pq._division_base(
        _W["div_units"][d], _W["flow_by_div"].get(d), _W["cap_km"], d)
    alpha = _W["alpha"]
    d_arr = alpha * fterm + (1.0 - alpha) * gterm

    Z = linkage(d_arr, method="average")
    M = squareform(d_arr)

    total_n, k_list = _W["total_n"], _W["k_list"]
    kdiv_of, ref_of = {}, {}
    for k_nat in k_list:
        k_div = max(2, min(int(round(n / (total_n / k_nat))), n))
        kdiv_of[k_nat] = k_div
        ref_of[k_nat] = fcluster(Z, t=k_div, criterion="maxclust")

    rng = np.random.default_rng(_W["seed"] + int(d))
    m = max(10, int(round(_W["subfrac"] * n)))
    acc = defaultdict(lambda: defaultdict(list))   # k -> ref cluster -> jaccards

    for _ in range(_W["nboot"]):
        S = np.sort(rng.choice(n, size=m, replace=False))
        Zs = linkage(squareform(M[np.ix_(S, S)], checks=False), method="average")
        for k_nat in k_list:
            k_div = min(kdiv_of[k_nat], m)
            if k_div < 2:
                continue
            rep = fcluster(Zs, t=k_div, criterion="maxclust")
            ids, jac = _max_jaccard(ref_of[k_nat][S], rep)
            for cid, jv in zip(ids, jac):
                acc[k_nat][cid].append(jv)

    out = []
    for k_nat in k_list:
        per = acc.get(k_nat)
        if not per:
            continue
        means = np.array([np.mean(v) for v in per.values()])
        sizes = np.array([np.sum(ref_of[k_nat] == c) for c in per.keys()],
                         dtype=float)
        out.append((k_nat, n, kdiv_of[k_nat], len(means),
                    float(means.mean()),
                    float((means * sizes).sum() / sizes.sum()),
                    float((means >= 0.75).mean()),
                    float((means <= 0.50).mean())))
    return d, n, out


def main(level, divlookup, pflow, cap_km, alpha, seed, out_csv):
    kgrid = os.environ.get("HS_KGRID", "")
    if not kgrid:
        sys.exit("hc_stability: HS_KGRID not set")
    k_list = [float(x) for x in kgrid.split(",")]
    nboot = int(os.environ.get("HS_NBOOT", "50"))
    subfrac = float(os.environ.get("HS_SUBFRAC", "0.75"))

    samp = os.environ.get("HS_SAMPLE_DIVS", "")
    sample_set = {int(x) for x in samp.replace(",", " ").split()} if samp else None

    div_units, code_of, dvA, dvB, total_n = pq._load_divisions(divlookup)
    print(f"  [{level}] {len(div_units)} divisions, {total_n:,} units; "
          f"k grid {k_list[0]:g}..{k_list[-1]:g} ({len(k_list)} values); "
          f"B={nboot} subfrac={subfrac} alpha={alpha} cap={cap_km}", flush=True)
    flow_by_div = pq._load_flows_by_div(pflow, code_of, dvA, dvB)

    work = [d for d in sorted(div_units)
            if sample_set is None or d in sample_set]
    nproc = max(1, min(int(os.environ.get("HS_NPROC", os.cpu_count() or 1)),
                       len(work)))
    print(f"  [{level}] bootstrapping {len(work)} divisions across "
          f"{nproc} workers ...", flush=True)

    initargs = (div_units, flow_by_div, cap_km, alpha, k_list, total_n,
                nboot, subfrac, seed)
    if nproc == 1:
        _init(*initargs)
        results = (_process_div(d) for d in work)
        pool = None
    else:
        pool = Pool(processes=nproc, initializer=_init, initargs=initargs)
        results = pool.imap_unordered(_process_div, work)

    rows = defaultdict(list)
    for d, n, div_rows in results:
        print(f"  div {d}: n={n}, cells={len(div_rows)}", flush=True)
        for r in div_rows:
            rows[r[0]].append(r)
    if pool is not None:
        pool.close(); pool.join()

    out = []
    for k_nat in sorted(rows):
        v = np.array([[r[1], r[3], r[4], r[5], r[6], r[7]]
                      for r in rows[k_nat]], dtype=float)
        w = v[:, 0]
        out.append({
            "geo_level": level, "alpha": alpha, "k_national": k_nat,
            "method": "HC",
            "jaccard_mean":  float((v[:, 2] * w).sum() / w.sum()),
            "jaccard_wmean": float((v[:, 3] * w).sum() / w.sum()),
            "frac_stable":   float((v[:, 4] * w).sum() / w.sum()),
            "frac_dissolved": float((v[:, 5] * w).sum() / w.sum()),
            "n_clusters": int(v[:, 1].sum()),
            "n_div_used": len(v), "n_units": int(w.sum()),
            "n_boot": nboot, "subfrac": subfrac,
        })
    pd.DataFrame(out).to_csv(out_csv, index=False)
    print(f"  wrote {len(out)} rows to {out_csv}", flush=True)


if __name__ == "__main__":
    if len(sys.argv) != 8:
        sys.exit(__doc__)
    main(level=sys.argv[1], divlookup=sys.argv[2], pflow=sys.argv[3],
         cap_km=float(sys.argv[4]), alpha=float(sys.argv[5]),
         seed=int(sys.argv[6]), out_csv=sys.argv[7])
