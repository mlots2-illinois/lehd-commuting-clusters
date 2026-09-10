#!/usr/bin/env bash
#===============================================================================
# fetch_lodes.sh — Reassemble the raw public inputs for the LEHD clustering
#                  pipeline (LODES OD flows + OMB CBSA delineation).
#
# All inputs fetched here are public and free from the U.S. Census Bureau. This
# script downloads the LEHD-LODES v8.4 Origin-Destination files for BOTH job-type
# universes by default -- JT00 ("all jobs", the baseline) and JT01 ("primary
# jobs", used by the job-type robustness rerun) -- for the 49 analyzed
# jurisdictions (lower-48 states + DC) and years 2019-2023, plus the 2023 OMB
# metropolitan delineation used for the CBSA preservation check. Set
# JOB_TYPES="JT00" to fetch the baseline universe only and halve the download.
#
# TIGER/Line shapefiles and 2020 Centers of Population are large and versioned;
# their source URLs are printed at the end for manual download into
# 00_shapefiles/ and 01_raw/ (see README.md at the repository root).
#
# Usage:
#     bash fetch_lodes.sh [DEST_DIR]
#
# DEST_DIR must be the pipeline's raw-input folder, i.e. <project>/03_data/01_raw,
# because the pipeline reads LODES from $raw_lodes = $proj/03_data/01_raw/lodes.
# If DEST_DIR is omitted, LEHD_PROJ is used when exported; otherwise the script
# stops and asks, rather than depositing ~70 GB somewhere the pipeline will not
# look for it.
#
# Downloaded layout (matches what data_02_create_files.do expects):
#     DEST_DIR/lodes/<st>_od_main_<JT>_<year>.csv.gz     JT in {JT00, JT01}
#     DEST_DIR/lodes/<st>_od_aux_<JT>_<year>.csv.gz
#     DEST_DIR/cbsa/list1_2023.xlsx
#
# Requires: curl (or wget), ~70 GB free disk (~35 GB for JT00 alone).
#===============================================================================
set -euo pipefail

# ---------------------------------------------------------------------------
# Configuration (matches _inc/baseline_namespaces.do and _master.do)
# ---------------------------------------------------------------------------
LODES_VER="LODES8"            # v8.x series on the LEHD server
# Job-type universes to fetch. JT00 ("all jobs") is the baseline behavioral
# input used throughout the paper; JT01 ("primary jobs") is consumed by the
# build phase (build_01_pflows loops over both) and feeds the inline job-type
# robustness rerun (do_robust_jobtype in _master.do -> the jobtype_sweep table).
# To skip JT01 entirely (halves the download; baseline analysis unaffected), set
#   JOB_TYPES="JT00"
# and drop JT01 from `job_types` in _master.do (set beside `levels`).
JOB_TYPES="${JOB_TYPES:-JT00 JT01}"
YEAR_START=2019
YEAR_END=2023
BASE="https://lehd.ces.census.gov/data/lodes/${LODES_VER}"
OMB_URL="https://www2.census.gov/programs-surveys/metro-micro/geographies/reference-files/2023/delineation-files/list1_2023.xlsx"

# Destination: explicit argument wins; otherwise derive it from LEHD_PROJ. We
# deliberately do NOT fall back to a path beside this script -- that is not where
# the pipeline reads LODES from, and silently downloading ~70 GB to the wrong
# place is far more costly than stopping here.
if [[ $# -ge 1 ]]; then
    DEST="$1"
elif [[ -n "${LEHD_PROJ:-}" ]]; then
    DEST="${LEHD_PROJ}/03_data/01_raw"
    echo "Using LEHD_PROJ: destination = ${DEST}"
else
    cat >&2 <<'EOF'
ERROR: no destination given and LEHD_PROJ is not exported.

The pipeline reads LODES from  <project>/03_data/01_raw/lodes , so pass that
raw-input folder explicitly:

    bash fetch_lodes.sh /path/to/project/03_data/01_raw

or export the project root once and re-run with no argument:

    export LEHD_PROJ=/path/to/project
    bash fetch_lodes.sh

See README.md at the repository root for the full input layout.
EOF
    exit 2
fi

LODES_DIR="${DEST}/lodes"
CBSA_DIR="${DEST}/cbsa"
mkdir -p "$LODES_DIR" "$CBSA_DIR"

# 49 analyzed jurisdictions: lower-48 states + DC (drops AK 02, HI 15, and
# territories 60/66/69/72/78, exactly as build_01_pflows.do does).
STATES=(al az ar ca co ct de dc fl ga id il in ia ks ky la me md ma mi mn ms \
        mo mt ne nv nh nj nm ny nc nd oh ok or pa ri sc sd tn tx ut vt va wa \
        wv wi wy)

# ---------------------------------------------------------------------------
# Downloader: prefer curl, fall back to wget; skip files already present.
# ---------------------------------------------------------------------------
dl () {  # dl <url> <outfile>
    local url="$1" out="$2"
    if [[ -s "$out" ]]; then
        echo "  skip (exists): $(basename "$out")"; return 0
    fi
    if command -v curl >/dev/null 2>&1; then
        curl -fsSL --retry 3 -o "$out" "$url" || { echo "  MISSING on server: $url" >&2; rm -f "$out"; return 1; }
    else
        wget -q -t 3 -O "$out" "$url" || { echo "  MISSING on server: $url" >&2; rm -f "$out"; return 1; }
    fi
    echo "  got: $(basename "$out")"
}

# ---------------------------------------------------------------------------
# 1. LODES Origin-Destination files (main + aux), per state per year
# ---------------------------------------------------------------------------
echo "Downloading LODES ${LODES_VER} OD files [${JOB_TYPES}], ${YEAR_START}-${YEAR_END}, to ${LODES_DIR}"
fail=0
for jt in $JOB_TYPES; do
    for st in "${STATES[@]}"; do
        for yr in $(seq "$YEAR_START" "$YEAR_END"); do
            for part in main aux; do
                f="${st}_od_${part}_${jt}_${yr}.csv.gz"
                dl "${BASE}/${st}/od/${f}" "${LODES_DIR}/${f}" || fail=$((fail+1))
            done
        done
    done
done
echo "LODES download complete (${fail} file(s) unavailable — some state-years may be absent upstream)."

# ---------------------------------------------------------------------------
# 2. OMB 2023 metropolitan delineation (for the CBSA crosswalk, data_04)
# ---------------------------------------------------------------------------
echo "Downloading OMB 2023 CBSA delineation to ${CBSA_DIR}"
dl "$OMB_URL" "${CBSA_DIR}/list1_2023.xlsx" || echo "  fetch list1_2023.xlsx manually (see README.md)."

# ---------------------------------------------------------------------------
# 3. Large/versioned inputs — pointers only (download manually)
# ---------------------------------------------------------------------------
cat <<'EOF'

-------------------------------------------------------------------------------
Manual downloads still required (large / versioned; see README.md):

  * 2020 TIGER/Line shapefiles (tract + county + state)  ->  00_shapefiles/
      https://www2.census.gov/geo/tiger/TIGER2020/TRACT/
      https://www2.census.gov/geo/tiger/TIGER2020/COUNTY/
      https://www2.census.gov/geo/tiger/TIGER2020/STATE/
    Unpack to  00_shapefiles/us_{tract,county,state}_2020/US_{...}_2020.shp

  * 2020 Centers of Population, TRACT file only   ->  01_raw/
      https://www.census.gov/geographies/reference-files/time-series/geo/
        centers-population.2020.html
    Download CenPop2020_Mean_TR.txt and RENAME it to
      01_raw/tract_populationcenters_2020.csv
    (already comma-delimited; renaming is the only step). The county centers
    file is NOT needed -- county centroids are derived from the tract centres.

  * County -> commuting-zone crosswalk (USDA/ERS; optional)  ->  01_raw/czones/
      county_to_cz.dta or .csv, columns: geoid (5-digit county FIPS), cz
    Needed only by assess_06_czones.do (tab:cz_comparison). Not a Census
    product; the pipeline skips that one table cleanly if it is absent.

After fetching, set $proj, $program and the Python path at the top of
code/_master.do (or export LEHD_PROJ / LEHD_PROG / LEHD_PY), then run the
pipeline. See README.md.
-------------------------------------------------------------------------------
EOF
echo "Done."
