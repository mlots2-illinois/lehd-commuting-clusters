"""
fetch_acs_pairs.py -- downloads the wider set of tract-level ACS 5-year aggregates (twenty counts from twelve tables) plus Gazetteer land area for the multi-relationship unit-of-analysis work.

Called from assess/assess_11_robinson.do as:
    python fetch_acs_pairs.py <out_csv>
    also accepted: --year <yyyy> (default 2023) and --force
Reads:  https://www2.census.gov/programs-surveys/acs/summary_file/<year>/table-based-SF/data/5YRData/acsdt5y<year>-<table>.dat for the tables in SPEC;
        https://www2.census.gov/geo/docs/maps-data/data/gazetteer/<year>_Gazetteer/<year>_Gaz_tracts_national.zip (all downloaded at run time)
Writes: <out_csv> with geoid, the SPEC columns (sorted) and aland (square metres); the raw files under <dir of out_csv>/acs_sf_cache/, reused when present
Notes:  companion to fetch_acs.py, which pulls the four tables of the single relationship. Every quantity is an ACS aggregate or count, never a median or a published rate.
        That is deliberate: the analysis aggregates tracts into arbitrary areal units, and only sums may be aggregated; rates and means are formed after aggregation from summed numerators and denominators.
        A construct whose table line is missing for a tract is written as missing rather than short-summed. Line numbers were verified against published national 2019-2023 totals; see assess_11_robinson.do.
        Needs network access at run time; the caller runs it only when <out_csv> is missing.
"""

import sys
import os
import time
import csv
import argparse
import urllib.request
import zipfile

BASE = ("https://www2.census.gov/programs-surveys/acs/summary_file/"
        "{year}/table-based-SF/data/5YRData")
GAZ = ("https://www2.census.gov/geo/docs/maps-data/data/gazetteer/"
       "{year}_Gazetteer/{year}_Gaz_tracts_national.zip")
TRACT_PREFIX = "1400000US"

# output column -> (ACS table, [lines to sum])
#
# Line numbers were verified against published national 2019-2023 totals
# before use; see the header comment of assess_11_robinson.do for the checks.
SPEC = {
    "pop":       ("b05002", ["B05002_E001"]),   # total population
    "foreign":   ("b05002", ["B05002_E013"]),   # foreign born
    "ed_tot":    ("b15003", ["B15003_E001"]),   # population 25+
    "ed_lths":   ("b15003", ["B15003_E%03d" % i for i in range(2, 17)]),
    "ed_ba":     ("b15003", ["B15003_E%03d" % i for i in range(22, 26)]),
    "workers":   ("b08301", ["B08301_E001"]),   # workers 16+
    "transit":   ("b08301", ["B08301_E010"]),   # public transportation
    "pov_univ":  ("b17001", ["B17001_E001"]),   # poverty universe
    "pov_below": ("b17001", ["B17001_E002"]),   # below poverty
    "hu":        ("b25002", ["B25002_E001"]),   # housing units
    "vacant":    ("b25002", ["B25002_E003"]),   # vacant units
    "occ":       ("b25003", ["B25003_E001"]),   # occupied units
    "renter":    ("b25003", ["B25003_E003"]),   # renter-occupied
    "agg_rent":  ("b25065", ["B25065_E001"]),   # aggregate gross rent
    "clf":       ("b23025", ["B23025_E003"]),   # civilian labor force
    "unemp":     ("b23025", ["B23025_E005"]),   # unemployed
    "agg_hhinc": ("b19025", ["B19025_E001"]),   # aggregate household income
    "hh":        ("b11001", ["B11001_E001"]),   # households
    "agg_earn":  ("b19313", ["B19313_E001"]),   # aggregate earnings
    "n_earn":    ("b20001", ["B20001_E001"]),   # population with earnings
}


def download(url, dest, retries=4):
    last = None
    for attempt in range(retries):
        try:
            urllib.request.urlretrieve(url, dest)
            if os.path.getsize(dest) > 0:
                return
            last = "empty file"
        except Exception as exc:
            last = exc
        time.sleep(2.0 * (attempt + 1))
    raise RuntimeError("download failed for {}: {}".format(url, last))


def fetch_table(year, table, cachedir, cols):
    name = "acsdt5y{}-{}.dat".format(year, table)
    dest = os.path.join(cachedir, name)
    if not os.path.exists(dest) or os.path.getsize(dest) == 0:
        download("{}/{}".format(BASE.format(year=year), name), dest)
    out = {}
    with open(dest, "r", encoding="utf-8") as fh:
        header = fh.readline().rstrip("\n").split("|")
        idx = {n: i for i, n in enumerate(header)}
        missing = [c for c in cols if c not in idx]
        if missing:
            raise RuntimeError("{}: columns absent: {}".format(table, missing))
        geo_i = idx["GEO_ID"]
        want = [(c, idx[c]) for c in cols]
        for raw in fh:
            rec = raw.rstrip("\n").split("|")
            gid = rec[geo_i]
            if not gid.startswith(TRACT_PREFIX):
                continue
            row = out.setdefault(gid[len(TRACT_PREFIX):], {})
            for c, j in want:
                txt = rec[j].strip()
                if not txt:
                    continue
                try:
                    val = float(txt)
                except ValueError:
                    continue
                if val >= 0:          # negative codes are jam values
                    row[c] = val
    print("  {}: {} tract rows".format(table, len(out)))
    return out


def fetch_land(year, cachedir):
    """Tract land area (ALAND, square metres) from the national gazetteer."""
    zpath = os.path.join(cachedir, "gaz_tracts_{}.zip".format(year))
    if not os.path.exists(zpath) or os.path.getsize(zpath) == 0:
        download(GAZ.format(year=year), zpath)
    land = {}
    with zipfile.ZipFile(zpath) as zf:
        member = [n for n in zf.namelist() if n.lower().endswith(".txt")][0]
        with zf.open(member) as fh:
            head = fh.readline().decode("latin-1").rstrip("\n").split("\t")
            head = [h.strip() for h in head]
            gi, ai = head.index("GEOID"), head.index("ALAND")
            for raw in fh:
                rec = raw.decode("latin-1").rstrip("\n").split("\t")
                try:
                    land[rec[gi].strip()] = float(rec[ai])
                except (ValueError, IndexError):
                    continue
    print("  gazetteer: {} tract rows".format(len(land)))
    return land


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("out_csv")
    ap.add_argument("--year", default="2023")
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    if os.path.exists(args.out_csv) and not args.force:
        print("fetch_acs_pairs: {} already exists; skipping (use --force)."
              .format(args.out_csv))
        return

    cachedir = os.path.join(os.path.dirname(os.path.abspath(args.out_csv)),
                            "acs_sf_cache")
    if not os.path.isdir(cachedir):
        os.makedirs(cachedir)

    tables = sorted({tab for tab, _ in SPEC.values()})
    data = {}
    for tab in tables:
        cols = sorted({c for _k, (t, cs) in SPEC.items() if t == tab
                       for c in cs})
        for geoid, vals in fetch_table(args.year, tab, cachedir, cols).items():
            data.setdefault(geoid, {}).update(vals)

    land = fetch_land(args.year, cachedir)

    out_cols = sorted(SPEC)
    tmp = args.out_csv + ".tmp"
    with open(tmp, "w", newline="", encoding="utf-8") as fh:
        wr = csv.writer(fh)
        wr.writerow(["geoid"] + out_cols + ["aland"])
        for geoid in sorted(data):
            row = data[geoid]
            cells = []
            for key in out_cols:
                lines = SPEC[key][1]
                vals = [row.get(c) for c in lines]
                # a table line missing for this tract makes the whole
                # construct missing rather than silently short-summed
                cells.append("" if any(v is None for v in vals)
                             else repr(sum(vals)))
            cells.append(repr(land[geoid]) if geoid in land else "")
            wr.writerow([geoid] + cells)
    os.replace(tmp, args.out_csv)
    print("fetch_acs_pairs: wrote {} tract rows -> {}"
          .format(len(data), args.out_csv))


if __name__ == "__main__":
    main()
