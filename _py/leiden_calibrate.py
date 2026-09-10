"""
leiden_calibrate.py -- Leiden resolution sweep on the sampled divisions: per (division, resolution) cell, cluster count, modularity, cross-seed stability (mean pairwise NMI), largest-community and singleton fractions.

Called from calib/calib_02_resolution.do as:
    LEIDEN_NSEEDS=<n> python leiden_calibrate.py <dissim_dir> <out_csv> <res1,res2,...> <seed> cpm
    (calib_02 passes $leiden_nseeds_<lvl> = 10; the fifth argument selects CPM or the RB configuration model)
Reads:  <dissim_dir>/dissim_div*.csv (the sampled divisions only)
Writes: <out_csv>, one row per (div_id, resolution)
Notes:  each cell runs N_SEEDS Leiden partitions (seed + 0 ... N-1) and reports means over seeds; stability is the mean NMI over all seed pairs, so N = 10 gives 45 comparisons per cell.
        The graph for a division is cached per worker process across resolutions. LEIDEN_NPROC caps the pool.
"""

import os
import sys
import time
from itertools import combinations
from multiprocessing import Pool
from pathlib import Path
from statistics import fmean

import leidenalg as la
import pandas as pd
from sklearn.metrics import normalized_mutual_info_score

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _lehd_common as common

DIV_RE = common.DIV_RE
N_SEEDS = int(os.environ.get("LEIDEN_NSEEDS", "5"))

_GRAPH_CACHE: dict = {}

def _partition_class(partition_type: str):
    if partition_type == "rb":
        return la.RBConfigurationVertexPartition
    return la.CPMVertexPartition

def _get_graph(csv_path):
    key = str(csv_path)
    cached = _GRAPH_CACHE.get(key)
    if cached is None:
        cached = common.build_similarity_graph(csv_path)
        _GRAPH_CACHE[key] = cached
    return cached

def _process_cell(args):
    csv_path, res, seed, partition_type = args
    pcls = _partition_class(partition_type)
    div_id, _verts, g = _get_graph(csv_path)
    nv = g.vcount()

    memberships = []
    mods = []
    ncl = []
    lfrac = []
    sfrac = []
    for s in range(N_SEEDS):
        part = la.find_partition(
            g, pcls, weights="weight",
            resolution_parameter=res, seed=seed + s,
        )
        memberships.append(list(part.membership))
        mods.append(part.modularity)
        ncl.append(len(part))
        lfrac.append(max(len(c) for c in part) / nv if nv else 1.0)
        n_singletons = sum(1 for c in part if len(c) == 1)
        sfrac.append(n_singletons / nv if nv else 1.0)

    if len(memberships) > 1:
        stability = fmean(
            normalized_mutual_info_score(a, b)
            for a, b in combinations(memberships, 2)
        )
    else:
        stability = 1.0

    row = {
        "div_id":       div_id,
        "resolution":   res,
        "n_vertices":   nv,
        "n_edges":      g.ecount(),
        "n_clusters":   fmean(ncl),
        "modularity":   fmean(mods),
        "stability":    stability,
        "largest_frac": fmean(lfrac),
        "singleton_frac": fmean(sfrac),
    }
    log_line = (
        f"  div {div_id} res={res:g}: "
        f"n_clusters={fmean(ncl):.1f} "
        f"stability={stability:.3f} "
        f"largest_frac={fmean(lfrac):.3f} "
        f"singleton_frac={fmean(sfrac):.3f}"
    )
    return row, log_line

def main(dissim_dir: str, out_csv: str, resolutions: list[float], seed: int,
         partition_type: str = "cpm") -> None:
    in_dir = Path(dissim_dir)
    csv_files = sorted(in_dir.glob("dissim_div*.csv"))
    if not csv_files:
        sys.exit(f"No dissim_div*.csv files found in {in_dir}")

    tasks = [(cp, res, seed, partition_type)
             for cp in csv_files for res in resolutions]

    nproc = int(os.environ.get("LEIDEN_NPROC", os.cpu_count() or 1))
    nproc = max(1, min(nproc, len(tasks)))
    print(f"Calibrating Leiden on {len(csv_files)} division(s) x "
          f"{len(resolutions)} resolutions = {len(tasks)} cells across "
          f"{nproc} workers ({N_SEEDS} seeds each) ...", flush=True)

    total = len(tasks)
    start = time.time()

    def _progress(done, line):
        el = time.time() - start
        print(f"  [{done:>{len(str(total))}}/{total} {100*done/total:4.0f}% "
              f"{el:6.1f}s] {line.strip()}",
              flush=True)

    rows: list[dict] = []
    done = 0
    if nproc == 1:
        for t in tasks:
            row, line = _process_cell(t)
            done += 1
            _progress(done, line)
            rows.append(row)
    else:
        with Pool(processes=nproc) as pool:
            for row, line in pool.imap_unordered(_process_cell, tasks):
                done += 1
                _progress(done, line)
                rows.append(row)
    print(f"  ... {total} cells done in {time.time() - start:.1f}s", flush=True)

    df = pd.DataFrame(rows).sort_values(
        ["div_id", "resolution"], kind="stable").reset_index(drop=True)
    df.to_csv(out_csv, index=False)
    print(f"Wrote {len(df)} rows to {out_csv} (partition={partition_type})",
          flush=True)

if __name__ == "__main__":
    if len(sys.argv) not in (5, 6):
        sys.exit("usage: leiden_calibrate.py <dissim_dir> <out_csv> "
                 "<res1,res2,...> <seed> [cpm|rb]")
    dissim_dir = sys.argv[1]
    out_csv    = sys.argv[2]
    resolutions = [float(x) for x in sys.argv[3].split(",")]
    seed = int(sys.argv[4])
    ptype = sys.argv[5] if len(sys.argv) == 6 else "cpm"
    main(dissim_dir, out_csv, resolutions, seed, ptype)
