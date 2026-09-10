"""
hc_run.py -- average-linkage hierarchical clustering (HC) of each division's dissimilarity matrix, cut at a per-division k, written as one assignment table.

Called from cluster/cluster_01_hc.do and assess/assess_01_elbow.do as:
    python hc_run.py <dissim_dir> <out_csv> --k-map "div:k,div:k,..."
and from calib/calib_03_hc_k.do as:
    python hc_run.py <dissim_dir> <out_csv> --k-list "div:k1|k2|...,div:..."
Reads:  <dissim_dir>/dissim_div*.csv
Writes: <out_csv> with div_id, geoid, cut_value (the requested k), hc_cluster
Notes:  cuts are by cluster count (--k-map or --k-list); --cut-heights, or a positional list, cuts by dendrogram height instead.
        k is capped at the division's unit count. Units whose D is 1 to every other unit are dropped before linkage (_lehd_common.load_dissim_arrays_connected). HC_NPROC caps the worker pool.
"""
from __future__ import annotations

import argparse
import os
import sys
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _lehd_common as common

DIV_RE = common.DIV_RE

def _build_linkage(csv_path: Path) -> tuple[int, list[str], np.ndarray]:
    div_id, tracts, i_arr, j_arr, d_arr = common.load_dissim_arrays_connected(csv_path)
    n = len(tracts)

    assert (i_arr < j_arr).all(), \
        f"div {div_id}: upper-triangle contract violated (need geo_i < geo_j)"
    # The condensed vector is indexed by pair rank, so it relies on the file holding every pair; a missing pair would read as D = 0 (identical), not as unrelated.
    cond = np.zeros(n * (n - 1) // 2)
    cond[i_arr * (2 * n - i_arr - 1) // 2 + (j_arr - i_arr - 1)] = d_arr

    Z = linkage(cond, method="average")
    return div_id, tracts, Z

def process_division_heights(csv_path: Path, cut_heights: list[float]) -> list[pd.DataFrame]:
    div_id, tracts, Z = _build_linkage(csv_path)
    n = len(tracts)
    frames = []
    counts: list[int] = []
    for h in cut_heights:
        labels = fcluster(Z, t=h, criterion="distance")
        counts.append(len(set(labels)))
        frames.append(pd.DataFrame({
            "div_id":     div_id,
            "geoid":      tracts,
            "cut_value":  h,
            "hc_cluster": labels.astype(int),
        }))
    print(f"  div {div_id:3d}: n={n}, "
          f"clusters at h={cut_heights[0]}->{cut_heights[-1]}: {counts}")
    return frames

def process_division_k(csv_path: Path, k: int) -> list[pd.DataFrame]:
    div_id, tracts, Z = _build_linkage(csv_path)
    n = len(tracts)
    # A division cannot hold more clusters than units; fcluster may also return fewer than k when tied merge heights make the requested count unattainable.
    k_eff = min(k, n)
    labels = fcluster(Z, t=k_eff, criterion="maxclust")
    actual_k = len(set(labels))
    print(f"  div {div_id:3d}: n={n}, requested k={k}, effective k={k_eff}, actual k={actual_k}")
    return [pd.DataFrame({
        "div_id":     div_id,
        "geoid":      tracts,
        "cut_value":  float(k),
        "hc_cluster": labels.astype(int),
    })]

def process_division_k_list(csv_path: Path, k_list: list[int]) -> list[pd.DataFrame]:
    div_id, tracts, Z = _build_linkage(csv_path)
    n = len(tracts)
    frames: list[pd.DataFrame] = []
    seen: set[int] = set()
    counts: list[tuple[int, int]] = []
    for k in k_list:
        # calib_03 passes a k grid that can repeat after proportional allocation to small divisions; each distinct k is cut once.
        if k in seen:
            continue
        seen.add(k)
        k_eff = min(k, n)
        labels = fcluster(Z, t=k_eff, criterion="maxclust")
        counts.append((k, len(set(labels))))
        frames.append(pd.DataFrame({
            "div_id":     div_id,
            "geoid":      tracts,
            "cut_value":  float(k),
            "hc_cluster": labels.astype(int),
        }))
    print(f"  div {div_id:3d}: n={n}, k->actual: {counts}")
    return frames

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_dir")
    parser.add_argument("out_csv")
    parser.add_argument("legacy_cuts", nargs="?", default=None,
                        help="Legacy positional: comma-separated cut heights")
    parser.add_argument("--cut-heights", default=None,
                        help="Comma-separated dendrogram heights")
    parser.add_argument("--k-map", default=None,
                        help='Per-division k as "div:k,div:k,..."')
    parser.add_argument("--k-list", default=None,
                        help='Per-division k list as "div:k1|k2|...,div:k1|k2|..."')
    args = parser.parse_args()

    in_dir = Path(args.input_dir)
    csv_files = sorted(in_dir.glob("dissim_div*.csv"))
    if not csv_files:
        sys.exit(f"No dissim_div*.csv files found in {in_dir}")

    cut_heights = None
    k_map: dict[int, int] | None = None
    k_list_map: dict[int, list[int]] | None = None
    if args.k_list is not None:
        k_list_map = {}
        for tok in args.k_list.split(","):
            d, ks = tok.split(":")
            k_list_map[int(d)] = [int(x) for x in ks.split("|")]
    elif args.k_map is not None:
        k_map = {}
        for tok in args.k_map.split(","):
            d, k = tok.split(":")
            k_map[int(d)] = int(k)
    elif args.cut_heights is not None:
        cut_heights = [float(x) for x in args.cut_heights.split(",")]
    # Height-based cuts, for use outside the pipeline.
    elif args.legacy_cuts is not None:
        cut_heights = [float(x) for x in args.legacy_cuts.split(",")]
    else:
        sys.exit("Must specify --cut-heights, --k-map, --k-list, or legacy positional cut heights")

    tasks: list[tuple] = []
    for cp in csv_files:
        if k_list_map is not None:
            div_id = common.parse_div_id(cp)
            if div_id not in k_list_map:
                print(f"  div {div_id}: no k in --k-list, skipping")
                continue
            tasks.append(("kl", cp, k_list_map[div_id]))
        elif k_map is not None:
            div_id = common.parse_div_id(cp)
            if div_id not in k_map:
                print(f"  div {div_id}: no k in --k-map, skipping")
                continue
            tasks.append(("k", cp, k_map[div_id]))
        else:
            tasks.append(("h", cp, cut_heights))

    nproc = int(os.environ.get("HC_NPROC", os.cpu_count() or 1))
    nproc = max(1, min(nproc, len(tasks)))
    print(f"Processing {len(tasks)} divisions across {nproc} workers ...")

    frames: list[pd.DataFrame] = []
    if nproc == 1:
        for t in tasks:
            frames.extend(_run_task(t))
    else:
        with Pool(processes=nproc) as pool:
            for result in pool.imap_unordered(_run_task, tasks):
                frames.extend(result)

    out = pd.concat(frames, ignore_index=True)
    out.to_csv(args.out_csv, index=False)
    print(f"Wrote {len(out)} rows to {args.out_csv}")

def _run_task(task):
    mode, cp, payload = task
    if mode == "kl":
        return process_division_k_list(cp, payload)
    if mode == "k":
        return process_division_k(cp, payload)
    return process_division_heights(cp, payload)

if __name__ == "__main__":
    main()
