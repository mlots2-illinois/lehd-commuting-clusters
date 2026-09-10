"""
partition_quality.py -- modularity, silhouette and mean within-cluster D of a given partition, scored per division on the stored dissimilarity files and averaged with unit weights.

Called from assess/assess_07_partition_quality.do as:
    python partition_quality.py <dissim_dir> <assign_csv> <out_csv>
Reads:  <dissim_dir>/dissim_div*.csv ($hc_input/<lvl>_<run_tag>, left by cluster_01_hc.do); <assign_csv> with columns geoid, cluster
Writes: <out_csv>, one row: modularity, silhouette, within_d, n_div_used, n_units
Notes:  unlike partition_quality_alpha.py this scores any labeling (HC or Leiden, any run) but needs the dissim files on disk. Every unit in a division must carry a label or the script exits.
        Divisions with fewer than 3 units, one cluster, or all singletons return nothing and are left out of the weighted mean. The silhouette needs the full n x n matrix per division. PQ_NPROC caps the pool.
"""

import os
import sys
from pathlib import Path
from multiprocessing import Pool

import numpy as np
import pandas as pd
from sklearn.metrics import silhouette_score

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _lehd_common as common

DIV_RE = common.DIV_RE

def _eval_division(csv_path: Path, label_of: dict):
    _, verts, i_arr, j_arr, d_arr = common.load_dissim_arrays(csv_path)
    n = len(verts)
    if n < 3:
        return None

    missing = [g for g in verts if g not in label_of]
    if missing:
        sys.exit(f"{len(missing)} geoids in {csv_path.name} have no cluster "
                 f"label, e.g. {missing[:5]}")

    raw = [label_of[g] for g in verts]
    codes, _ = pd.factorize(pd.Series(raw))
    k = len(set(codes))
    if k < 2 or k >= n:
        return None

    M = np.ones((n, n), dtype=np.float64)
    np.fill_diagonal(M, 0.0)
    M[i_arr, j_arr] = d_arr
    M[j_arr, i_arr] = d_arr

    sil = float(silhouette_score(M, codes, metric="precomputed"))

    g = common.similarity_graph(n, i_arr, j_arr, d_arr)
    mod = float(g.modularity(codes.tolist(), weights="weight"))

    codes_arr = np.asarray(codes)
    iu_i, iu_j = np.triu_indices(n, k=1)
    same = codes_arr[iu_i] == codes_arr[iu_j]
    within_d = float(M[iu_i[same], iu_j[same]].mean()) if same.any() else float("nan")

    return n, mod, sil, within_d

_LABEL_OF = None

def _init(label_of):
    global _LABEL_OF
    _LABEL_OF = label_of

def _eval_worker(csv_path):
    return _eval_division(csv_path, _LABEL_OF)

def main(dissim_dir, assign_csv, out_csv):
    a = pd.read_csv(assign_csv, dtype={"geoid": str})
    a.columns = [c.strip().lower() for c in a.columns]
    label_of = dict(zip(a["geoid"], a["cluster"]))

    files = [f for f in sorted(Path(dissim_dir).glob("dissim_div*.csv"))
             if DIV_RE.search(f.name)]

    nproc = int(os.environ.get("PQ_NPROC", os.cpu_count() or 1))
    nproc = max(1, min(nproc, len(files))) if files else 1
    print(f"Scoring {len(files)} divisions across {nproc} workers ...", flush=True)

    rows = []
    if nproc == 1:
        _init(label_of)
        results = (_eval_worker(f) for f in files)
    else:
        pool = Pool(processes=nproc, initializer=_init, initargs=(label_of,))
        results = pool.imap_unordered(_eval_worker, files)
    for res in results:
        if res is not None:
            rows.append(res)
    if nproc != 1:
        pool.close()
        pool.join()

    if not rows:
        pd.DataFrame([{"modularity": float("nan"), "silhouette": float("nan"),
                       "within_d": float("nan"), "n_div_used": 0,
                       "n_units": 0}]).to_csv(out_csv, index=False)
        return

    arr = np.array(rows, dtype=float)
    nvec = arr[:, 0]
    tot = nvec.sum()
    out = {
        "modularity":  float((arr[:, 1] * nvec).sum() / tot),
        "silhouette":  float((arr[:, 2] * nvec).sum() / tot),
        "within_d":    float((arr[:, 3] * nvec).sum() / tot),
        "n_div_used":  int(len(rows)),
        "n_units":     int(tot),
    }
    pd.DataFrame([out]).to_csv(out_csv, index=False)

if __name__ == "__main__":
    if len(sys.argv) != 4:
        sys.exit("usage: partition_quality.py <dissim_dir> <assign_csv> <out_csv>")
    main(sys.argv[1], sys.argv[2], sys.argv[3])
