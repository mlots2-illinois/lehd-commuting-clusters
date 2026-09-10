"""
calib_flow_coverage.py -- P-weighted share of flow pairs within each candidate distance cap, pooled and as a per-division 10/50/90 percentile band.

Called from calib/calib_01_alpha_cap.do as:
    python calib_flow_coverage.py <calib_flow.dta> <cap1,...> <out_csv>
Reads:  <calib_flow.dta> (div_id d_km P_ij), pairs with P_ij > 0 only
Writes: <out_csv> with cap_km, coverage, cov_p10, cov_p50, cov_p90
Notes:  coverage is the cumulative sum of P_ij over pairs sorted by distance, so it is monotone in the cap; the selector reads it as the share of commuting the cap does not truncate.
"""

import sys
import numpy as np
import pandas as pd

def main(flow_dta, caps, out_csv):
    df = pd.read_stata(flow_dta, columns=["div_id", "d_km", "P_ij"])
    df = df[df["P_ij"] > 0]
    caps = np.asarray(sorted(caps), dtype=float)

    d_all = df["d_km"].to_numpy(float)
    p_all = df["P_ij"].to_numpy(float)
    order = np.argsort(d_all, kind="stable")
    d_sorted = d_all[order]
    p_cum = np.cumsum(p_all[order])
    tot = p_cum[-1] if len(p_cum) else 0.0
    idx = np.searchsorted(d_sorted, caps, side="right")
    pooled = np.where(idx > 0, p_cum[np.clip(idx - 1, 0, len(p_cum) - 1)], 0.0)
    pooled = pooled / tot if tot > 0 else np.zeros_like(caps)

    per_div = []
    for _, sub in df.groupby("div_id"):
        d = sub["d_km"].to_numpy(float)
        p = sub["P_ij"].to_numpy(float)
        o = np.argsort(d, kind="stable")
        ds = d[o]
        pc = np.cumsum(p[o])
        t = pc[-1]
        if t <= 0:
            continue
        ii = np.searchsorted(ds, caps, side="right")
        cov = np.where(ii > 0, pc[np.clip(ii - 1, 0, len(pc) - 1)], 0.0) / t
        per_div.append(cov)
    if per_div:
        band = np.vstack(per_div)
        p10 = np.nanpercentile(band, 10, axis=0)
        p50 = np.nanpercentile(band, 50, axis=0)
        p90 = np.nanpercentile(band, 90, axis=0)
    else:
        p10 = p50 = p90 = np.full_like(caps, np.nan)

    out = pd.DataFrame({
        "cap_km": caps,
        "coverage": pooled,
        "cov_p10": p10,
        "cov_p50": p50,
        "cov_p90": p90,
    })
    out.to_csv(out_csv, index=False)
    print(f"  wrote {len(out)} caps to {out_csv} "
          f"(pooled coverage {pooled.min():.3f}..{pooled.max():.3f})", flush=True)

if __name__ == "__main__":
    if len(sys.argv) != 4:
        sys.exit("usage: calib_flow_coverage.py <calib_flow_dta> <cap_csv> <out_csv>")
    cap_list = [float(x) for x in sys.argv[2].replace(",", " ").split()]
    main(sys.argv[1], cap_list, sys.argv[3])
