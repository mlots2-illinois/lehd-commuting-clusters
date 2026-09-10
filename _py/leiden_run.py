"""
leiden_run.py -- runs Leiden (CPM by default) on every division's similarity graph at each resolution, best of N seeds, and writes the memberships.

Called from cluster/cluster_02_leiden.do as:
    python leiden_run.py <dissim_dir> <out_csv> <res1,res2,...> <seed> cpm
    (<res1,...> is cluster_02's res_grid; the fifth argument selects CPM or the RB configuration model)
Reads:  <dissim_dir>/dissim_div*.csv (one file per division)
Writes: <out_csv> with div_id, geoid, resolution, leiden_cluster
Notes:  N = LEIDEN_NSEEDS (default 5). cluster_02_leiden.do does not set it; calib_02_resolution.do sets 10 for leiden_calibrate.py.
        Seeds are <seed> + 0 ... N-1 and the partition with the highest quality function wins. LEIDEN_NPROC caps the worker pool.
"""

import os
import sys
from multiprocessing import Pool
from pathlib import Path

import leidenalg as la
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _lehd_common as common

DIV_RE = common.DIV_RE
N_SEEDS = int(os.environ.get("LEIDEN_NSEEDS", "5"))

def _partition_class(partition_type: str):
    if partition_type == "rb":
        return la.RBConfigurationVertexPartition
    return la.CPMVertexPartition

def _build_graph(csv_path: Path):
    return common.build_similarity_graph(csv_path)

# Leiden is stochastic; we keep the best of N_SEEDS runs on quality() (the CPM or RB objective), not on modularity, so the choice matches what the partition optimises.
def _best_partition(g, pcls, res, base_seed):
    best = None
    best_q = None
    for s in range(N_SEEDS):
        part = la.find_partition(
            g, pcls, weights="weight",
            resolution_parameter=res, seed=base_seed + s,
        )
        q = part.quality()
        if best_q is None or q > best_q:
            best, best_q = part, q
    return best

def _process_div(args):
    csv_path, resolutions, seed, partition_type = args
    pcls = _partition_class(partition_type)
    div_id, verts, g = _build_graph(csv_path)

    rows: list[pd.DataFrame] = []
    log_lines: list[str] = []
    for res in resolutions:
        part = _best_partition(g, pcls, res, seed)
        rows.append(pd.DataFrame({
            "div_id":         div_id,
            "geoid":          verts,
            "resolution":     res,
            "leiden_cluster": part.membership,
        }))
        log_lines.append(
            f"  div {div_id} res={res:g}: {len(part)} clusters, "
            f"modularity={part.modularity:.4f} (best of {N_SEEDS} seeds)"
        )
    return pd.concat(rows, ignore_index=True), log_lines

def main(dissim_dir: str, out_csv: str, resolutions: list[float], seed: int,
         partition_type: str = "cpm") -> None:
    in_dir = Path(dissim_dir)
    csv_files = sorted(in_dir.glob("dissim_div*.csv"))
    if not csv_files:
        sys.exit(f"No dissim_div*.csv files found in {in_dir}")

    tasks = [(cp, resolutions, seed, partition_type) for cp in csv_files]

    # One task per division; imap_unordered returns divisions as they finish, so the row order varies between runs while the memberships do not.
    nproc = int(os.environ.get("LEIDEN_NPROC", os.cpu_count() or 1))
    nproc = max(1, min(nproc, len(tasks)))
    print(f"Running Leiden on {len(tasks)} divisions across {nproc} workers "
          f"({N_SEEDS} seeds each) ...")

    frames: list[pd.DataFrame] = []
    if nproc == 1:
        for t in tasks:
            df, log = _process_div(t)
            for line in log:
                print(line)
            frames.append(df)
    else:
        with Pool(processes=nproc) as pool:
            for df, log in pool.imap_unordered(_process_div, tasks):
                for line in log:
                    print(line)
                frames.append(df)

    out = pd.concat(frames, ignore_index=True)
    out.to_csv(out_csv, index=False)
    print(f"Wrote {len(out)} rows to {out_csv}")

if __name__ == "__main__":
    if len(sys.argv) not in (5, 6):
        sys.exit("usage: leiden_run.py <dissim_dir> <out_csv> "
                 "<res1,res2,...> <seed> [cpm|rb]")
    main(
        dissim_dir=sys.argv[1],
        out_csv=sys.argv[2],
        resolutions=[float(x) for x in sys.argv[3].split(",")],
        seed=int(sys.argv[4]),
        partition_type=sys.argv[5] if len(sys.argv) == 6 else "cpm",
    )
