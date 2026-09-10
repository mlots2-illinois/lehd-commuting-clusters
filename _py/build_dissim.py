"""
build_dissim.py -- writes one dissimilarity file per division (D for every unit pair) from the division lookup and the pooled flow file; the Python twin of _inc/build_div_dissim.do.

Called from cluster/cluster_01_hc.do as:
    python build_dissim.py <divlookup_dta> <pflow_dta> <outdir> <alpha> <cap_km> <divfmt> <counts_csv>
    options: --divs d1,d2,... to restrict to a subset of divisions; --nproc N for the worker pool
Reads:  <divlookup_dta> (div_id geoid lat lon); <pflow_dta> (geo_i geo_j P_ij)
Writes: <outdir>/dissim_div<NNNNN>.csv with geo_i, geo_j, D (stale files deleted first); <counts_csv> with div_id, n_units
Notes:  D = alpha * (1 - min(P_ij, 1)) + (1 - alpha) * min(d_km / cap_km, 1), with haversine distance in km; a pair with no flow row gets P_ij = 0.
        Only flow rows with geo_i < geo_j are read, matching the upper-triangle pair layout. DISSIM_NPROC caps the worker pool.
"""
from __future__ import annotations

import argparse
import math
import os
import re
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd

R_KM = 6371.0
_TORAD = math.pi / 180.0

def _divfmt_width(divfmt: str) -> int:
    m = re.search(r"%0?(\d+)", divfmt)
    return int(m.group(1)) if m else 5

def _build_one(task):
    (div_id, geoids, lat, lon, f_i, f_j, f_P, alpha, cap_km, out_path) = task
    n = len(geoids)
    if n < 2:
        return (div_id, n, False)

    iu, ju = np.triu_indices(n, k=1)

    la_i = lat[iu] * _TORAD
    la_j = lat[ju] * _TORAD
    dlat = (lat[ju] - lat[iu]) * _TORAD
    dlon = (lon[ju] - lon[iu]) * _TORAD
    a = np.sin(dlat / 2.0) ** 2 + np.cos(la_i) * np.cos(la_j) * np.sin(dlon / 2.0) ** 2
    d_km = 2.0 * R_KM * np.arctan2(np.sqrt(a), np.sqrt(1.0 - a))

    gi = geoids[iu]
    gj = geoids[ju]

    if len(f_i):
        flows = pd.DataFrame({"geo_i": f_i, "geo_j": f_j, "P_ij": f_P})
        pairs = pd.DataFrame({"geo_i": gi, "geo_j": gj})
        merged = pairs.merge(flows, on=["geo_i", "geo_j"], how="left")
        P = merged["P_ij"].fillna(0.0).to_numpy(dtype=np.float64)
    else:
        P = np.zeros(len(iu), dtype=np.float64)

    P_clip = np.minimum(P, 1.0)
    D = alpha * (1.0 - P_clip) + (1.0 - alpha) * np.minimum(d_km / cap_km, 1.0)

    out = pd.DataFrame({"geo_i": gi, "geo_j": gj, "D": D})
    out.to_csv(out_path, index=False)
    return (div_id, n, True)

def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("divlookup")
    ap.add_argument("pflow")
    ap.add_argument("outdir")
    ap.add_argument("alpha", type=float)
    ap.add_argument("cap_km", type=float)
    ap.add_argument("divfmt")
    ap.add_argument("counts_out")
    ap.add_argument("--divs", default=None,
                    help="Comma-separated division ids to restrict to (default: all)")
    ap.add_argument("--nproc", type=int, default=None)
    args = ap.parse_args()

    width = _divfmt_width(args.divfmt)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    for f in outdir.glob("dissim_div*.csv"):
        f.unlink()

    dl = pd.read_stata(args.divlookup, columns=["div_id", "geoid", "lat", "lon"])
    dl["geoid"] = dl["geoid"].astype(str)
    dl["div_id"] = dl["div_id"].astype(int)

    restrict = None
    if args.divs:
        restrict = {int(x) for x in args.divs.split(",")}
        dl = dl[dl["div_id"].isin(restrict)]

    pf = pd.read_stata(args.pflow, columns=["geo_i", "geo_j", "P_ij"])
    pf["geo_i"] = pf["geo_i"].astype(str)
    pf["geo_j"] = pf["geo_j"].astype(str)
    pf = pf[pf["geo_i"] < pf["geo_j"]]

    tasks = []
    for div_id, g in dl.groupby("div_id", sort=True):
        g = g.sort_values("geoid")
        geoids = g["geoid"].to_numpy()
        if len(geoids) < 2:
            print(f"  build_dissim div {div_id}: only {len(geoids)} unit(s), skipping")
            continue
        gset = set(geoids.tolist())
        sub = pf[pf["geo_i"].isin(gset) & pf["geo_j"].isin(gset)]
        fname = f"dissim_div{int(div_id):0{width}d}.csv"
        tasks.append((
            int(div_id),
            geoids,
            g["lat"].to_numpy(dtype=np.float64),
            g["lon"].to_numpy(dtype=np.float64),
            sub["geo_i"].to_numpy(),
            sub["geo_j"].to_numpy(),
            sub["P_ij"].to_numpy(dtype=np.float64),
            args.alpha,
            args.cap_km,
            str(outdir / fname),
        ))

    if not tasks:
        raise SystemExit("build_dissim: no divisions with >=2 units to build.")

    nproc = args.nproc or int(os.environ.get("DISSIM_NPROC", os.cpu_count() or 1))
    nproc = max(1, min(nproc, len(tasks)))
    print(f"Building {len(tasks)} divisions across {nproc} workers "
          f"(alpha={args.alpha}, cap_km={args.cap_km}) ...")

    counts = []
    if nproc == 1:
        for t in tasks:
            div_id, n, built = _build_one(t)
            if built:
                counts.append((div_id, n))
                print(f"  div {div_id}: {n} units -> {n*(n-1)//2:,} pairs")
    else:
        with Pool(processes=nproc) as pool:
            for div_id, n, built in pool.imap_unordered(_build_one, tasks):
                if built:
                    counts.append((div_id, n))
                    print(f"  div {div_id}: {n} units -> {n*(n-1)//2:,} pairs")

    counts_df = pd.DataFrame(sorted(counts), columns=["div_id", "n_units"])
    counts_df.to_csv(args.counts_out, index=False)
    print(f"Wrote {len(counts_df)} division CSVs; counts -> {args.counts_out}")

if __name__ == "__main__":
    main()
