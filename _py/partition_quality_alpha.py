"""
partition_quality_alpha.py -- self-containment, silhouette and modularity of the HC partition across an alpha grid (and optionally a k grid), rebuilt in memory from the division lookup and the pooled flows.

Called from calib/calib_01_alpha_cap.do (Phase D), calib/calib_03_hc_k.do and assess/assess_08_partition_quality_alpha.do as:
    python partition_quality_alpha.py <level> <divlookup_dta> <pflow_dta> <cap_km> <hc_k_national> <leiden_res> <seed> <alpha1,...> <out_csv>
    env: PQ_KGRID (k grid, replaces <hc_k_national>), PQ_SAMPLE_DIVS, PQ_LIGHT=1, PQ_FULL=1, PQ_MINFLOW, PQ_NPROC, PFLOW_CHUNK
Reads:  <divlookup_dta> (div_id geoid lat lon); <pflow_dta> (geo_i geo_j P_ij f_ij), read in chunks
Writes: <out_csv>, one row per (alpha, k_national); the k_national column is dropped when only one k is scored
Notes:  only the HC arm is scored and the scoring is deterministic; <leiden_res> and <seed> keep the call signature aligned with the calibration scripts.
        Light mode (PQ_LIGHT=1, or PQ_KGRID without PQ_FULL=1) skips the n x n silhouette matrix and returns modularity and silhouette as missing.
        hc_stability.py and null_structure_test.py import the loaders and _division_base from this module.
"""

import os
import sys
from collections import defaultdict
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage
from sklearn.metrics import silhouette_score

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _lehd_common as common

EARTH_KM = 6371.0
CHUNK = int(os.environ.get("PFLOW_CHUNK", "5000000"))

def _norm_geoid(s: pd.Series, width: int) -> pd.Series:
    s = s.astype(str).str.strip()
    s = s.str.replace(r"\.0$", "", regex=True)
    return s.str.zfill(width)

def _haversine(lat_i, lon_i, lat_j, lon_j):
    r = np.pi / 180.0
    a = (np.sin((lat_j - lat_i) * r / 2.0) ** 2
         + np.cos(lat_i * r) * np.cos(lat_j * r)
         * np.sin((lon_j - lon_i) * r / 2.0) ** 2)
    return 2.0 * EARTH_KM * np.arctan2(np.sqrt(a), np.sqrt(1.0 - a))

def _load_divisions(divlookup_dta):
    dl = pd.read_stata(divlookup_dta, columns=["div_id", "geoid", "lat", "lon"])
    dl["geoid"] = dl["geoid"].astype(str).str.strip()
    dl["geoid"] = dl["geoid"].str.replace(r"\.0$", "", regex=True)
    width = int(dl["geoid"].str.len().max())
    dl["geoid"] = dl["geoid"].str.zfill(width)
    dl["div_id"] = dl["div_id"].astype(int)
    dl = dl.drop_duplicates(["geoid", "div_id"])

    geoids = np.sort(dl["geoid"].unique())
    code_of = {g: i for i, g in enumerate(geoids)}
    dl["code"] = dl["geoid"].map(code_of).to_numpy()

    nunits = len(geoids)
    # A unit can sit in two overlapping divisions (core of one, margin of the next); dvA and dvB hold both so a flow is routed to every division containing both ends.
    dvA = np.full(nunits, -1, dtype=np.int32)
    dvB = np.full(nunits, -1, dtype=np.int32)
    latlon = dl.drop_duplicates("code").set_index("code")[["lat", "lon"]]

    div_units = {}
    for d, sub in dl.groupby("div_id"):
        codes = np.sort(sub["code"].to_numpy())
        for c in codes:
            if dvA[c] == -1:
                dvA[c] = d
            elif dvB[c] == -1 and dvA[c] != d:
                dvB[c] = d
        ll = latlon.loc[codes]
        div_units[d] = pd.DataFrame({"code": codes,
                                     "lat": ll["lat"].to_numpy(),
                                     "lon": ll["lon"].to_numpy()})
    return div_units, code_of, dvA, dvB, nunits

def _load_flows_by_div(pflow_dta, code_of, dvA, dvB):
    """Route within-division flows to their division.

    PQ_MINFLOW applies the same f_ij floor that build_01_pflows applies when
    it writes the _mf<m> pool files, so a minimum-flow robustness check can
    be run off the baseline pool without rebuilding it.
    """
    min_flow = float(os.environ.get("PQ_MINFLOW", "0"))
    acc = defaultdict(lambda: ([], [], [], []))

    def _route(ci, cj, P, f, da, db):
        for d in np.unique(da):
            if d < 0:
                continue
            m = da == d
            store = acc[d]
            store[0].append(ci[m]); store[1].append(cj[m])
            store[2].append(P[m]);  store[3].append(f[m])

    width = max(len(g) for g in code_of)

    reader = pd.read_stata(pflow_dta, columns=["geo_i", "geo_j", "P_ij", "f_ij"],
                           chunksize=CHUNK)
    nrows = 0
    for chunk in reader:
        nrows += len(chunk)
        ci = _norm_geoid(chunk["geo_i"], width).map(code_of).to_numpy()
        cj = _norm_geoid(chunk["geo_j"], width).map(code_of).to_numpy()
        valid = ~(pd.isna(ci) | pd.isna(cj))
        match = valid.mean()
        # A geoid dtype mismatch (lost leading zeros) zeroes the flow term, so we stop here.
        if match < 0.99:
            sys.exit(f"only {match:.1%} of pflow geoids matched divlookup; "
                     f"dtype/zero-padding mismatch?")
        if not valid.any():
            continue
        ci = ci[valid].astype(np.int64); cj = cj[valid].astype(np.int64)
        P = chunk["P_ij"].to_numpy(float)[valid]
        f = chunk["f_ij"].to_numpy(float)[valid]
        if min_flow > 0:
            keep = f >= min_flow
            ci, cj, P, f = ci[keep], cj[keep], P[keep], f[keep]
            if not keep.any():
                continue
        iA, iB = dvA[ci], dvB[ci]
        jA, jB = dvA[cj], dvB[cj]
        for da, db in ((iA, jA), (iA, jB), (iB, jA), (iB, jB)):
            m = (da == db) & (da >= 0)
            if m.any():
                _route(ci[m], cj[m], P[m], f[m], da[m], db[m])

    out = {}
    for d, (li, lj, lp, lf) in acc.items():
        out[d] = (np.concatenate(li), np.concatenate(lj),
                  np.concatenate(lp), np.concatenate(lf))
    print(f"  [flows] scanned {nrows:,} pflow rows (min_flow={min_flow:g}); "
          f"{len(out)} divisions have within-division flow", flush=True)
    return out

def _triu_rank(i, j, n):
    return i * (2 * n - i - 1) // 2 + (j - i - 1)

def _division_base(sub, flows, cap_km, div_id):
    codes = sub["code"].to_numpy()
    n = len(codes)
    lat = sub["lat"].to_numpy(float)
    lon = sub["lon"].to_numpy(float)
    iu, ju = np.triu_indices(n, k=1)
    gterm = np.minimum(_haversine(lat[iu], lon[iu], lat[ju], lon[ju]) / cap_km, 1.0)
    fterm = np.ones(len(iu), dtype=np.float64)
    raw_flow = np.zeros(len(iu), dtype=np.float64)
    raw_prop = np.zeros(len(iu), dtype=np.float64)

    if flows is not None:
        ci, cj, P, f = flows
        li = np.searchsorted(codes, ci)
        lj = np.searchsorted(codes, cj)
        a = np.minimum(li, lj)
        b = np.maximum(li, lj)
        rank = _triu_rank(a, b, n)
        assert len(np.unique(rank)) == len(rank), \
            f"div {div_id}: duplicate flow pairs"
        fterm[rank] = 1.0 - np.minimum(P, 1.0)
        raw_flow[rank] = f
        raw_prop[rank] = P
    return n, iu, ju, fterm, gterm, raw_flow, raw_prop

def _score(codes, M, i_idx, j_idx, d_arr):
    n = M.shape[0]
    k = len(set(codes))
    if k < 2 or k >= n:
        return None
    sil = float(silhouette_score(M, codes, metric="precomputed"))
    g = common.similarity_graph(n, i_idx, j_idx, d_arr)
    mod = float(g.modularity(list(codes), weights="weight"))
    codes_arr = np.asarray(codes)
    same = codes_arr[i_idx] == codes_arr[j_idx]
    within_d = float(d_arr[same].mean()) if same.any() else float("nan")
    return mod, sil, within_d

def _self_containment(codes, i_idx, j_idx, w):
    """Within-cluster share of a pair weight.

    Called with raw f_ij for the volume-weighted score and with P_ij for the
    proportional-flow score.  P_ij is the quantity the dissimilarity is
    actually built on, so the proportional score is the one that measures
    what the clustering optimises; the raw score is reported alongside it.
    """
    same = codes[i_idx] == codes[j_idx]
    return float(w[same].sum()), float(w.sum())

_W = {}

def _init(div_units, flow_by_div, cap_km, alphas, k_list, total_n, light):
    _W.update(div_units=div_units, flow_by_div=flow_by_div, cap_km=cap_km,
              alphas=alphas, k_list=k_list, total_n=total_n, light=light)

def _process_div(d):
    div_units, flow_by_div = _W["div_units"], _W["flow_by_div"]
    cap_km, alphas, k_list = _W["cap_km"], _W["alphas"], _W["k_list"]
    total_n, light = _W["total_n"], _W["light"]

    n = len(div_units[d])
    if n < 3:
        return d, n, []
    n, i_idx, j_idx, fterm, gterm, raw_flow, raw_prop = \
        _division_base(div_units[d], flow_by_div.get(d), cap_km, d)

    out = []
    for alpha in alphas:
        d_arr = alpha * fterm + (1.0 - alpha) * gterm
        Z = linkage(d_arr, method="average")
        if light:
            M = None
        else:
            M = np.ones((n, n), dtype=np.float64)
            np.fill_diagonal(M, 0.0)
            M[i_idx, j_idx] = d_arr
            M[j_idx, i_idx] = d_arr

        for k_nat in k_list:
            # The national k target is allocated to divisions in proportion to their unit count, the rule cluster_01_hc.do uses, so the scores describe the partition actually built.
            k_div = max(2, min(int(round(n / (total_n / k_nat))), n))
            hc = fcluster(Z, t=k_div, criterion="maxclust")
            if len(set(hc)) < 2 or len(set(hc)) >= n:
                continue
            num, den = _self_containment(hc, i_idx, j_idx, raw_flow)
            pnum, pden = _self_containment(hc, i_idx, j_idx, raw_prop)
            same = hc[i_idx] == hc[j_idx]
            within_d = float(d_arr[same].mean()) if same.any() else float("nan")
            if light:
                mod = sil = float("nan")
            else:
                s = _score(hc, M, i_idx, j_idx, d_arr)
                mod, sil = (s[0], s[1]) if s is not None else (float("nan"),) * 2
            out.append(((alpha, k_nat),
                        (n, mod, sil, within_d, len(set(hc)), num, den,
                         pnum, pden)))
    return d, n, out

def main(level, divlookup, pflow, cap_km, hc_k_national, leiden_res, seed,
         alphas, out_csv):
    kgrid_env = os.environ.get("PQ_KGRID", "")
    k_list = [float(x) for x in kgrid_env.split(",")] if kgrid_env else [hc_k_national]
    light = ((bool(kgrid_env) or os.environ.get("PQ_LIGHT", "0") == "1")
             and os.environ.get("PQ_FULL", "0") != "1")

    samp_env = os.environ.get("PQ_SAMPLE_DIVS", "")
    sample_set = {int(x) for x in samp_env.replace(",", " ").split()} if samp_env else None

    div_units, code_of, dvA, dvB, total_n = _load_divisions(divlookup)
    if sample_set is not None:
        kept = sorted(d for d in div_units if d in sample_set)
        print(f"  [{level}] restricting to {len(kept)} sampled divisions: {kept}",
              flush=True)
    print(f"  [{level}] {len(div_units)} divisions, {total_n:,} distinct units; "
          f"cap={cap_km} k_grid={k_list} light={light} alphas={alphas}", flush=True)
    flow_by_div = _load_flows_by_div(pflow, code_of, dvA, dvB)

    rows = defaultdict(list)

    work = [d for d in sorted(div_units)
            if sample_set is None or d in sample_set]
    nproc = int(os.environ.get("PQ_NPROC", os.cpu_count() or 1))
    nproc = max(1, min(nproc, len(work))) if work else 1
    print(f"  [{level}] scoring {len(work)} divisions across {nproc} workers ...",
          flush=True)

    initargs = (div_units, flow_by_div, cap_km, alphas, k_list, total_n, light)
    if nproc == 1:
        _init(*initargs)
        results = (_process_div(d) for d in work)
    else:
        pool = Pool(processes=nproc, initializer=_init, initargs=initargs)
        results = pool.imap_unordered(_process_div, work)
    for d, n, div_rows in results:
        print(f"  div {d}: n={n}, cells={len(div_rows)}", flush=True)
        for key, row in div_rows:
            rows[key].append(row)
    if nproc != 1:
        pool.close()
        pool.join()

    # Modularity, silhouette and within_d are unit-weighted means over divisions; self-containment instead sums numerators and denominators so it is a true flow share.
    def _wmean(vals_arr, wts):
        m = ~np.isnan(vals_arr)
        if not m.any():
            return float("nan"), 0
        return float((vals_arr[m] * wts[m]).sum() / wts[m].sum()), int(m.sum())

    out = []
    for (alpha, k_nat), vals in sorted(rows.items()):
        arr = np.array([(v[0], v[1], v[2], v[3]) for v in vals], dtype=float)
        nvec = arr[:, 0]
        tot = nvec.sum()
        sc_num = sum(v[5] for v in vals)
        sc_den = sum(v[6] for v in vals)
        sp_num = sum(v[7] for v in vals)
        sp_den = sum(v[8] for v in vals)
        mod_w, mod_nd = _wmean(arr[:, 1], nvec)
        sil_w, sil_nd = _wmean(arr[:, 2], nvec)
        wd_w,  wd_nd  = _wmean(arr[:, 3], nvec)
        out.append({
            "geo_level": level,
            "alpha": alpha,
            "k_national": k_nat,
            "method": "HC",
            "self_containment_prop": (sp_num / sp_den) if sp_den > 0 else float("nan"),
            "self_containment": (sc_num / sc_den) if sc_den > 0 else float("nan"),
            "modularity": mod_w,
            "silhouette": sil_w,
            "within_d": wd_w,
            "n_clusters": int(sum(v[4] for v in vals)),
            "n_div_used": len(vals),
            "n_units": int(tot),
            "n_div_modularity": mod_nd,
            "n_div_silhouette": sil_nd,
            "n_div_within_d": wd_nd,
        })
    df = pd.DataFrame(out)
    if len(k_list) == 1:
        df = df.drop(columns=["k_national"])
    df.to_csv(out_csv, index=False)
    print(f"  wrote {len(out)} rows to {out_csv}", flush=True)

if __name__ == "__main__":
    if len(sys.argv) != 10:
        sys.exit("usage: partition_quality_alpha.py <level> <divlookup_dta> "
                 "<pflow_dta> <cap_km> <hc_k_national> <leiden_res> <seed> "
                 "<alpha1,...> <out_csv>")
    main(level=sys.argv[1], divlookup=sys.argv[2], pflow=sys.argv[3],
         cap_km=float(sys.argv[4]),
         hc_k_national=float(sys.argv[5]), leiden_res=float(sys.argv[6]),
         seed=int(sys.argv[7]),
         alphas=[float(x) for x in sys.argv[8].split(",")],
         out_csv=sys.argv[9])
