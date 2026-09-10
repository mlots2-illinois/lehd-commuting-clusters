"""
_lehd_common.py -- shared readers that turn a dissim_div<NN>.csv file into index arrays or an igraph similarity graph.

Imported by hc_run.py, leiden_run.py, leiden_calibrate.py, partition_quality.py and partition_quality_alpha.py; not called from Stata.
Reads:  dissim_div<NN>.csv with columns geo_i, geo_j, D (geo_i < geo_j, every pair present)
Notes:  three edge-filter variants exist. load_dissim_arrays keeps every unit and every pair (partition_quality.py).
        load_dissim_arrays_connected drops units whose D is 1 to every other unit (hc_run.py).
        build_similarity_graph keeps the pairs with w = 1 - D > 0 and only the units on them (leiden_run.py, leiden_calibrate.py).
        similarity_graph builds the same edges on a caller-supplied vertex count, so isolated units keep their index (partition_quality.py, partition_quality_alpha.py).
"""
from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd
import igraph as ig

DIV_RE = re.compile(r"dissim_div(\d+)\.csv$")

def parse_div_id(csv_path: Path) -> int:
    m = DIV_RE.search(Path(csv_path).name)
    if m is None:
        raise ValueError(
            f"Filename does not match dissim_div<NN>.csv: {csv_path}")
    return int(m.group(1))

def read_dissim(csv_path: Path) -> tuple[int, pd.DataFrame]:
    div_id = parse_div_id(csv_path)
    df = pd.read_csv(csv_path, dtype={"geo_i": str, "geo_j": str},
                     usecols=["geo_i", "geo_j", "D"])
    return div_id, df

def load_dissim_arrays(csv_path: Path):
    div_id, df = read_dissim(csv_path)
    verts = sorted(set(df["geo_i"]).union(df["geo_j"]))
    idx = {g: i for i, g in enumerate(verts)}
    i_arr = df["geo_i"].map(idx).to_numpy(dtype=np.int64)
    j_arr = df["geo_j"].map(idx).to_numpy(dtype=np.int64)
    d_arr = df["D"].to_numpy(dtype=np.float64)
    return div_id, verts, i_arr, j_arr, d_arr

# HC drops units with no positive similarity so its vertex set matches the graph Leiden sees; the two methods then cluster the same units.
def load_dissim_arrays_connected(csv_path: Path):
    div_id, df = read_dissim(csv_path)
    w = 1.0 - df["D"].to_numpy(dtype=np.float64)
    connected = set(df.loc[w > 0.0, "geo_i"]).union(df.loc[w > 0.0, "geo_j"])
    keep = df["geo_i"].isin(connected) & df["geo_j"].isin(connected)
    df = df.loc[keep]
    verts = sorted(connected)
    idx = {g: i for i, g in enumerate(verts)}
    i_arr = df["geo_i"].map(idx).to_numpy(dtype=np.int64)
    j_arr = df["geo_j"].map(idx).to_numpy(dtype=np.int64)
    d_arr = df["D"].to_numpy(dtype=np.float64)
    return div_id, verts, i_arr, j_arr, d_arr

# Takes n from the caller so a unit with no positive-similarity edge stays a vertex; modularity is then computed over the same units as the silhouette.
def similarity_graph(n: int, i_arr: np.ndarray, j_arr: np.ndarray,
                     d_arr: np.ndarray) -> ig.Graph:
    w = 1.0 - np.asarray(d_arr, dtype=np.float64)
    keep = w > 0.0
    edges = np.column_stack((np.asarray(i_arr)[keep], np.asarray(j_arr)[keep]))
    g = ig.Graph(n=n, edges=edges, directed=False)
    g.es["weight"] = w[keep]
    return g

# w = 1 - D is the edge weight; a pair at the distance cap with no flow has D = 1 and carries no edge, which keeps the CPM graph sparse.
def build_similarity_graph(csv_path: Path):
    div_id, df = read_dissim(csv_path)
    w = 1.0 - df["D"].to_numpy(dtype=np.float64)
    keep = w > 0.0
    df = df.loc[keep]
    w = w[keep]

    verts = sorted(set(df["geo_i"]).union(df["geo_j"]))
    idx = {v: i for i, v in enumerate(verts)}
    edges = np.column_stack((df["geo_i"].map(idx).to_numpy(dtype=np.int64),
                             df["geo_j"].map(idx).to_numpy(dtype=np.int64)))
    g = ig.Graph(n=len(verts), edges=edges, directed=False)
    g.es["weight"] = w
    return div_id, verts, g
