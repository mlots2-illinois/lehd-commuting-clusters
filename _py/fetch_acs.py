"""
fetch_acs.py -- downloads the tract-level ACS 5-year aggregates for the earnings-on-education relationship and writes them as one CSV.

Called from assess/assess_09_maup.do and assess/assess_10_maup_paramsweep.do as:
    python fetch_acs.py <out_csv>
    also accepted: --year <yyyy> (default 2023) and --force
Reads:  https://www2.census.gov/programs-surveys/acs/summary_file/<year>/table-based-SF/data/5YRData/acsdt5y<year>-<table>.dat for b19313, b20001, b15003, b01003 (downloaded at run time)
Writes: <out_csv> with geoid, agg_earn, n_earn, ed_tot, ed_ba, pop; the raw .dat files under <dir of out_csv>/acs_sf_cache/
Notes:  needs network access on first run; the callers run it only when <out_csv> is missing, and a .dat file already in acs_sf_cache is reused.
        ed_ba sums bachelor's through doctorate (B15003 lines 22-25); negative ACS codes (jam values) become missing.
"""

import sys
import os
import time
import argparse
import urllib.request

BASE = ("https://www2.census.gov/programs-surveys/acs/summary_file/"
        "{year}/table-based-SF/data/5YRData")

NEEDED = {
    "b19313": ["B19313_E001"],
    "b20001": ["B20001_E001"],
    "b15003": ["B15003_E001", "B15003_E022", "B15003_E023",
               "B15003_E024", "B15003_E025"],
    "b01003": ["B01003_E001"],
}

TRACT_PREFIX = "1400000US"

def to_num(x):
    x = x.strip()
    if x == "":
        return None
    try:
        v = float(x)
    except ValueError:
        return None
    if v < 0:
        return None
    return v

def download(url, dest, retries=4):
    last = None
    for attempt in range(retries):
        try:
            urllib.request.urlretrieve(url, dest)
            size = os.path.getsize(dest)
            if size > 0:
                return size
            last = "empty file"
        except Exception as e:
            last = e
        time.sleep(2.0 * (attempt + 1))
    raise RuntimeError("download failed for {}: {}".format(url, last))

def fetch_table(year, table, cachedir):
    url = "{}/acsdt5y{}-{}.dat".format(BASE.format(year=year), year, table)
    dest = os.path.join(cachedir, "acsdt5y{}-{}.dat".format(year, table))
    if not (os.path.exists(dest) and os.path.getsize(dest) > 0):
        download(url, dest)
    cols = NEEDED[table]
    out = {}
    with open(dest, "r", encoding="utf-8") as f:
        header = f.readline().rstrip("\n").split("|")
        idx = {name: i for i, name in enumerate(header)}
        want = [(c, idx[c]) for c in cols]
        geo_i = idx["GEO_ID"]
        for raw in f:
            rec = raw.rstrip("\n").split("|")
            gid = rec[geo_i]
            if not gid.startswith(TRACT_PREFIX):
                continue
            geoid = gid[len(TRACT_PREFIX):]
            out[geoid] = {c: to_num(rec[j]) for c, j in want}
    print("  {}: {} tract rows".format(table, len(out)))
    return out

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("out_csv")
    ap.add_argument("--year", default="2023")
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    if os.path.exists(args.out_csv) and not args.force:
        print("fetch_acs: {} already exists; skipping (use --force to refetch)."
              .format(args.out_csv))
        return

    cachedir = os.path.join(os.path.dirname(os.path.abspath(args.out_csv)),
                            "acs_sf_cache")
    os.makedirs(cachedir, exist_ok=True)
    tabs = {t: fetch_table(args.year, t, cachedir) for t in NEEDED}

    geoids = set()
    for t in tabs.values():
        geoids.update(t.keys())

    def get(tab, col, geoid):
        rec = tabs[tab].get(geoid)
        return None if rec is None else rec.get(col)

    def fmt(v):
        return "" if v is None else repr(v)

    n = 0
    with open(args.out_csv, "w") as f:
        f.write("geoid,agg_earn,n_earn,ed_tot,ed_ba,pop\n")
        for geoid in sorted(geoids):
            agg_earn = get("b19313", "B19313_E001", geoid)
            n_earn = get("b20001", "B20001_E001", geoid)
            ed_tot = get("b15003", "B15003_E001", geoid)
            ba = get("b15003", "B15003_E022", geoid)
            ma = get("b15003", "B15003_E023", geoid)
            prof = get("b15003", "B15003_E024", geoid)
            doc = get("b15003", "B15003_E025", geoid)
            pop = get("b01003", "B01003_E001", geoid)
            ed_ba = None
            if None not in (ba, ma, prof, doc):
                ed_ba = ba + ma + prof + doc
            f.write(",".join([
                geoid, fmt(agg_earn), fmt(n_earn),
                fmt(ed_tot), fmt(ed_ba), fmt(pop),
            ]) + "\n")
            n += 1

    print("fetch_acs: wrote {} tract rows -> {}".format(n, args.out_csv))

if __name__ == "__main__":
    main()
