* ensure_county_div_lookup.do -- defines and runs ensure_county_div_lookup, which builds the single-division county lookup when it does not exist.
* Included by: cluster/cluster_01_hc.do, cluster_02_leiden.do, cluster_03_compile.do, cluster_05_reconcile.do, report/report_02_dissim_example.do.
* Expects: $temp
* Reads:  $temp/centroids_county.dta
* Writes: $temp/division_lookup_county.dta (only when missing)
* Notes:  counties are few enough to cluster nationally in one piece, so the county lookup is one division (div_id = 1, is_core = 1) and build_03_divisions.do is not needed at this level. The file is never rebuilt once present; delete it to refresh.

capture program drop ensure_county_div_lookup
program define ensure_county_div_lookup

    local outfile "$temp/division_lookup_county.dta"
    local src     "$temp/centroids_county.dta"

    local rebuild 0
    capture confirm file "`outfile'"
    if _rc local rebuild 1

    if `rebuild' {
        di as text "Building `outfile' from `src' ..."
        use geoid lat lon using "`src'", clear
        gen int  div_id  = 1
        gen byte is_core = 1
        order div_id geoid lat lon is_core
        sort div_id geoid

        label variable div_id   "Division ID (county-level: always 1)"
        label variable geoid    "5-digit county FIPS"
        label variable lat      "County centroid latitude"
        label variable lon      "County centroid longitude"
        label variable is_core  "1 = core of window (always 1 for county)"

        compress
        save "`outfile'", replace
    }
end

ensure_county_div_lookup
