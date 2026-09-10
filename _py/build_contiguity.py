"""
build_contiguity.py -- writes the Queen-contiguity edge list of a polygon shapefile.

Called from cluster/cluster_07_contiguity.do as:
    python build_contiguity.py <shp_file> <out_csv> GEOID
Reads:  <shp_file>
Writes: <out_csv> with geoid_i, geoid_j (each pair once, geoid_i < geoid_j)
Notes:  needs geopandas and libpysal (see requirements.txt). Islands (polygons with no neighbour) are counted and printed, not written.
"""
from __future__ import annotations

import sys

import geopandas as gpd
import pandas as pd
from libpysal.weights import Queen

def main(shp_path: str, out_csv: str, geoid_col: str = "GEOID") -> None:
    print(f"Reading shapefile: {shp_path}")
    gdf = gpd.read_file(shp_path)
    if geoid_col not in gdf.columns:
        if geoid_col.lower() in gdf.columns:
            geoid_col = geoid_col.lower()
        else:
            sys.exit(
                f"Column '{geoid_col}' not found.  Available: {list(gdf.columns)}"
            )

    n = len(gdf)
    print(f"  {n} polygons loaded")

    print("Computing Queen contiguity ...")
    w = Queen.from_dataframe(gdf, ids=gdf[geoid_col].tolist(), use_index=False)

    pairs = []
    for i, neighbors in w.neighbors.items():
        for j in neighbors:
            if str(i) < str(j):
                pairs.append((i, j))

    df = pd.DataFrame(pairs, columns=["geoid_i", "geoid_j"])
    df = df.drop_duplicates().sort_values(["geoid_i", "geoid_j"])
    df.to_csv(out_csv, index=False)
    print(f"Wrote {len(df)} edges to {out_csv}")

    n_islands = sum(1 for nbrs in w.neighbors.values() if len(nbrs) == 0)
    if n_islands:
        print(f"  Note: {n_islands} islands (no Queen neighbors)")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        sys.exit("usage: build_contiguity.py <shp_file> <out_csv> [<geoid_col>]")
    main(
        shp_path=sys.argv[1],
        out_csv=sys.argv[2],
        geoid_col=sys.argv[3] if len(sys.argv) > 3 else "GEOID",
    )
