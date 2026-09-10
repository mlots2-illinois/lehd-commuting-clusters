"""
random_partitions.py -- random areal partitions of the tract universe, matched on the number of units, as the null for assess_12_random_null.do.

Called from assess/assess_12_random_null.do as:
    python random_partitions.py --neighbors $temp/tract_neighbors.dta --universe $temp/final_assignments_<run_tag>_tract.dta --outdir <nulldir> --reps $null_reps --seed $seed
    also accepted: --k K1 K2 ... (default 644 3100 4298 6627)
Reads:  --universe (column geoid); --neighbors (two geoid columns, each edge once)
Writes: <outdir>/random_<unconstrained|contiguous>_k<K>.dta with geoid and one int32 column c1 ... c<reps>
Notes:  the unit-of-analysis result compares delineations built from commuting against the county map, and a sceptic can reply that any partition with a few thousand units would beat the county. Two nulls answer that.
        unconstrained assigns tracts to K clusters uniformly at random, so group means are near-unbiased draws from the tract distribution.
        contiguous grows regions by multi-source breadth-first search from K random seeds on the tract adjacency graph; groups are connected and compact like real geographies but blind to commuting, which is the reference that matters.
        Components without a seed become their own cluster, so the realised count can exceed K. The RNG is seeded with --seed + K (+1 for unconstrained), so each file is reproducible on its own.
"""

import os
import argparse
from collections import deque

import numpy as np
import pandas as pd


def load_universe(path):
    if path.endswith(".dta"):
        df = pd.read_stata(path, columns=["geoid"])
    else:
        df = pd.read_csv(path, dtype={"geoid": str}, usecols=["geoid"])
    return df["geoid"].astype(str).unique()


def build_adjacency(path, geoids):
    """Adjacency in CSR-like form, restricted to the analysis universe."""
    idx = {g: i for i, g in enumerate(geoids)}
    if path.endswith(".dta"):
        e = pd.read_stata(path)
    else:
        e = pd.read_csv(path, dtype=str)
    a = e.iloc[:, 0].astype(str).map(idx)
    b = e.iloc[:, 1].astype(str).map(idx)
    keep = a.notna() & b.notna()
    a = a[keep].to_numpy(np.int64)
    b = b[keep].to_numpy(np.int64)
    # symmetrise: the stored list may hold each pair once
    src = np.concatenate([a, b])
    dst = np.concatenate([b, a])
    order = np.argsort(src, kind="stable")
    src, dst = src[order], dst[order]
    n = len(geoids)
    start = np.searchsorted(src, np.arange(n + 1))
    return start, dst


def grow_contiguous(start, dst, n, k, rng):
    """Multi-source BFS from k random seeds; ties broken by queue order.

    Components containing no seed are handed their own cluster, so the
    realised cluster count can exceed k slightly on a disconnected graph.
    """
    lab = np.full(n, -1, np.int64)
    seeds = rng.choice(n, size=min(k, n), replace=False)
    q = deque()
    for c, s in enumerate(seeds):
        lab[s] = c
        q.append(s)
    while q:
        u = q.popleft()
        for v in dst[start[u]:start[u + 1]]:
            if lab[v] == -1:
                lab[v] = lab[u]
                q.append(v)
    # unreached tracts sit in components with no seed; give each component
    # its own cluster rather than dropping it
    nxt = int(lab.max()) + 1
    for i in np.flatnonzero(lab == -1):
        if lab[i] != -1:
            continue
        lab[i] = nxt
        q.append(i)
        while q:
            u = q.popleft()
            for v in dst[start[u]:start[u + 1]]:
                if lab[v] == -1:
                    lab[v] = nxt
                    q.append(v)
        nxt += 1
    return lab


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--neighbors", required=True)
    ap.add_argument("--universe", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--k", nargs="+", type=int,
                    default=[644, 3100, 4298, 6627])
    ap.add_argument("--reps", type=int, default=20)
    ap.add_argument("--seed", type=int, default=202600331)
    args = ap.parse_args()

    geoids = load_universe(args.universe)
    n = len(geoids)
    print("random_partitions: universe = %d tracts" % n)
    start, dst = build_adjacency(args.neighbors, geoids)
    deg = np.diff(start)
    print("  adjacency: %d directed edges, %d isolated tracts"
          % (len(dst), int((deg == 0).sum())))

    if not os.path.isdir(args.outdir):
        os.makedirs(args.outdir)

    for kind in ("unconstrained", "contiguous"):
        for k in args.k:
            rng = np.random.default_rng(args.seed + k +
                                        (0 if kind == "contiguous" else 1))
            out = {"geoid": geoids}
            realised = []
            for r in range(1, args.reps + 1):
                if kind == "unconstrained":
                    lab = rng.integers(0, k, size=n)
                else:
                    lab = grow_contiguous(start, dst, n, k, rng)
                out["c%d" % r] = lab.astype(np.int32)
                realised.append(len(np.unique(lab)))
            df = pd.DataFrame(out)
            path = os.path.join(args.outdir,
                                "random_%s_k%d.dta" % (kind, k))
            df.to_stata(path, write_index=False, version=118)
            print("  %-13s k=%-5d realised %d-%d clusters over %d reps -> %s"
                  % (kind, k, min(realised), max(realised), args.reps,
                     os.path.basename(path)))


if __name__ == "__main__":
    main()
