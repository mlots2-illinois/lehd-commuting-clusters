"""
calib_dispersion_grid.py -- dispersion statistics of D over the (alpha x cap) grid, per division, computed from the sampled pair table (calib_01 Phase C).

Called from calib/calib_01_alpha_cap.do as:
    python calib_dispersion_grid.py <calib_flow.dta> <alpha1,...> <cap1,...> <out_csv>
Reads:  <calib_flow.dta> (div_id d_km P_ij; the pairs table from _inc/build_div_dissim.do)
Writes: <out_csv> with alpha, cap_km, div_id, n_pairs, d_mean, d_sd, d_p25, d_p50, d_p75, d_min, d_max, d_frac_lo (D < 0.05), d_frac_hi (D > 0.95), iqr
Notes:  percentiles follow Stata's definition (_stata_pctiles) so the Python and Stata summaries agree to the digit. CALIB_NPROC caps the pool; each worker holds the whole pair table.
"""

import os
import sys
import time
from itertools import product
from multiprocessing import Pool

import numpy as np
import pandas as pd

_GROUPS = None

def _stata_pctiles(a_sorted, ps):
    n = a_sorted.size
    out = []
    for p in ps:
        q = n * p / 100.0
        i = int(np.floor(q))
        if q == i:
            if i <= 0:
                out.append(float(a_sorted[0]))
            elif i >= n:
                out.append(float(a_sorted[n - 1]))
            else:
                out.append(0.5 * (float(a_sorted[i - 1]) + float(a_sorted[i])))
        else:
            out.append(float(a_sorted[min(i, n - 1)]))
    return out

def _init(groups):
    global _GROUPS
    _GROUPS = groups

def _cell(args):
    alpha, cap = args
    rows = []
    for div_id, d_km, p_clip in _GROUPS:
        D = alpha * (1.0 - p_clip) + (1.0 - alpha) * np.minimum(d_km / cap, 1.0)
        n = D.size
        Ds = np.sort(D)
        p25, p50, p75 = _stata_pctiles(Ds, (25, 50, 75))
        rows.append({
            "alpha":     alpha,
            "cap_km":    cap,
            "div_id":    div_id,
            "n_pairs":   n,
            "d_mean":    float(D.mean()),
            "d_sd":      float(D.std(ddof=1)) if n > 1 else np.nan,
            "d_p25":     p25,
            "d_p50":     p50,
            "d_p75":     p75,
            "d_min":     float(Ds[0]),
            "d_max":     float(Ds[-1]),
            "d_frac_lo": float(np.mean(D < 0.05)),
            "d_frac_hi": float(np.mean(D > 0.95)),
            "iqr":       p75 - p25,
        })
    return rows

def _load_groups(flow_dta):
    df = pd.read_stata(flow_dta, columns=["div_id", "d_km", "P_ij"])
    df["div_id"] = df["div_id"].astype("int64")
    p_clip = np.minimum(df["P_ij"].to_numpy(dtype=np.float64), 1.0)
    d_km = df["d_km"].to_numpy(dtype=np.float64)
    div = df["div_id"].to_numpy()
    groups = []
    for d in np.unique(div):
        m = div == d
        groups.append((int(d), d_km[m], p_clip[m]))
    return groups

def main(flow_dta, alphas, caps, out_csv):
    groups = _load_groups(flow_dta)
    cells = list(product(alphas, caps))
    nproc = int(os.environ.get("CALIB_NPROC", os.cpu_count() or 1))
    nproc = max(1, min(nproc, len(cells)))
    total = len(cells)
    print(f"calib Phase C: {len(alphas)} alphas x {len(caps)} caps = "
          f"{total} cells over {len(groups)} division(s) "
          f"across {nproc} workers ...", flush=True)

    start = time.time()
    step = max(1, total // 20)

    def _beat(done):
        el = time.time() - start
        print(f"  [{done:>{len(str(total))}}/{total} {100*done/total:4.0f}% "
              f"{el:6.1f}s]", flush=True)

    rows = []
    done = 0
    if nproc == 1:
        _init(groups)
        for c in cells:
            rows.extend(_cell(c))
            done += 1
            if done % step == 0 or done == total:
                _beat(done)
    else:
        with Pool(processes=nproc, initializer=_init, initargs=(groups,)) as pool:
            for cell_rows in pool.imap_unordered(_cell, cells):
                rows.extend(cell_rows)
                done += 1
                if done % step == 0 or done == total:
                    _beat(done)

    cols = ["alpha", "cap_km", "div_id", "n_pairs", "d_mean", "d_sd",
            "d_p25", "d_p50", "d_p75", "d_min", "d_max",
            "d_frac_lo", "d_frac_hi", "iqr"]
    df = pd.DataFrame(rows, columns=cols).sort_values(
        ["alpha", "cap_km", "div_id"], kind="stable").reset_index(drop=True)
    df.to_csv(out_csv, index=False)
    print(f"Wrote {len(df)} rows to {out_csv} in {time.time() - start:.1f}s",
          flush=True)

if __name__ == "__main__":
    if len(sys.argv) != 5:
        sys.exit("usage: calib_dispersion_grid.py <calib_flow.dta> "
                 "<alpha1,...> <cap1,...> <out_sensitivity.csv>")
    flow_dta = sys.argv[1]
    alphas = [float(x) for x in sys.argv[2].split(",") if x != ""]
    caps = [float(x) for x in sys.argv[3].split(",") if x != ""]
    out_csv = sys.argv[4]
    main(flow_dta, alphas, caps, out_csv)
