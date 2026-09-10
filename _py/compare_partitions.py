"""
compare_partitions.py -- MI, NMI (arithmetic and geometric), AMI and ARI of two labelings, through scikit-learn.

Called from _inc/compare_partitions.do as:
    python compare_partitions.py <in_csv> <out_csv>
Reads:  <in_csv>, two columns of labels, one row per unit
Writes: <out_csv>, one row with mi, nmi_arith, nmi_geom, ami, ari, n
Notes:  labels are read as strings so numeric and string cluster ids compare alike. The Stata caller keeps ami, ari and nmi_geom; it computes its own MI, NMI and ARI in Stata.
"""

import sys

import pandas as pd
from sklearn.metrics import (
    adjusted_mutual_info_score,
    adjusted_rand_score,
    mutual_info_score,
    normalized_mutual_info_score,
)

def main(in_csv: str, out_csv: str) -> None:
    df = pd.read_csv(in_csv, dtype=str, keep_default_na=False)
    a = df.iloc[:, 0].to_numpy()
    b = df.iloc[:, 1].to_numpy()

    out = pd.DataFrame([{
        "mi":        mutual_info_score(a, b),
        "nmi_arith": normalized_mutual_info_score(a, b, average_method="arithmetic"),
        "nmi_geom":  normalized_mutual_info_score(a, b, average_method="geometric"),
        "ami":       adjusted_mutual_info_score(a, b, average_method="arithmetic"),
        "ari":       adjusted_rand_score(a, b),
        "n":         len(a),
    }])
    out.to_csv(out_csv, index=False)
    print(f"compare_partitions.py: n={len(a)} "
          f"ami={out['ami'].iloc[0]:.4f} ari={out['ari'].iloc[0]:.4f} "
          f"nmi_arith={out['nmi_arith'].iloc[0]:.4f}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit("usage: compare_partitions.py <in_csv> <out_csv>")
    main(sys.argv[1], sys.argv[2])
