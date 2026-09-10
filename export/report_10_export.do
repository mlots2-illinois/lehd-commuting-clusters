clear all
cap log close
local _logdate = subinstr("$S_DATE", " ", "", .)
log using "$log/report_10_export_${run_tag}_`_logdate'.log", replace text

* report_10_export.do -- export the final tract and county assignments with FIPS components, centroid coordinates and 2020 population as the distributable data product.
* Called by: _master.do (export phase, do_export switch).
* Reads:  $temp/final_assignments_<run_tag>_<lvl>.dta; $temp/centroids_<lvl>.dta; $temp/centroids_tract.dta
* Writes: $output/tract_clusters_<run_tag>.dta and .csv; $output/county_clusters_<run_tag>.dta and .csv
* Notes:  County population is summed from tract centroids so the two files share one population source. div_id and is_core are carried for reference; at the county level there is a single division.

foreach lvl of global levels {

    display as text _n "{hline 70}"
    display as text "===== report_10 LEVEL: `lvl' ====="
    display as text "{hline 70}"

    capture confirm file "$temp/final_assignments_${run_tag}_`lvl'.dta"
    if _rc {
        display as text "report_10: no final_assignments for `lvl'; skipping level."
        continue
    }
    capture confirm file "$temp/centroids_`lvl'.dta"
    if _rc {
        display as text "report_10: $temp/centroids_`lvl'.dta not found; skipping level."
        continue
    }
    capture confirm file "$temp/centroids_tract.dta"
    if _rc {
        display as text "report_10: $temp/centroids_tract.dta not found (needed for the pop column); skipping level `lvl'."
        continue
    }

    use "$temp/final_assignments_${run_tag}_`lvl'.dta", clear

    if "`lvl'" == "tract" {
        merge 1:1 geoid using "$temp/centroids_`lvl'.dta", ///
            keep(master match) keepusing(st_fips lat lon pop) nogenerate
        replace pop = 0 if mi(pop)
    }
    else {
        merge 1:1 geoid using "$temp/centroids_county.dta", ///
            keep(master match) keepusing(st_fips lat lon) nogenerate

        preserve
        use geoid pop using "$temp/centroids_tract.dta", clear
        gen str5 county_geoid = substr(geoid, 1, 5)
        gcollapse (sum) pop, by(county_geoid)
        rename county_geoid geoid
        tempfile cpops
        save `cpops'
        restore
        merge 1:1 geoid using `cpops', keep(master match) nogenerate
        replace pop = 0 if mi(pop)
    }

    if "`lvl'" == "tract" {
        gen str5 county_geoid = substr(geoid, 1, 5)
        gen str3 county_fips  = substr(geoid, 3, 3)

        label variable geoid           "11-digit tract FIPS"
        label variable st_fips         "State FIPS"
        label variable county_fips     "County FIPS (3-digit)"
        label variable county_geoid    "County FIPS (5-digit)"
        label variable lat             "Pop-weighted centroid latitude"
        label variable lon             "Pop-weighted centroid longitude"
        label variable pop             "2020 Census tract population"

        order geoid st_fips county_fips county_geoid lat lon pop ///
            hc_cluster leiden_cluster div_id is_core ///
            hc_was_fragment leiden_was_fragment

        local stem "tract_clusters_${run_tag}"
    }
    else {
        capture confirm variable st_fips
        if _rc gen str2 st_fips = substr(geoid, 1, 2)
        gen str3 county_fips = substr(geoid, 3, 3)

        label variable geoid           "5-digit county FIPS"
        label variable st_fips         "State FIPS"
        label variable county_fips     "County FIPS (3-digit)"
        label variable lat             "Centroid latitude"
        label variable lon             "Centroid longitude"
        label variable pop             "2020 Census county population"

        order geoid st_fips county_fips lat lon pop ///
            hc_cluster leiden_cluster div_id is_core ///
            hc_was_fragment leiden_was_fragment

        local stem "county_clusters_${run_tag}"
    }

    label variable hc_cluster      "Hierarchical clustering cluster ID"
    label variable leiden_cluster  "Leiden community detection cluster ID"
    label variable div_id          "Snake-scan division ID (county: always 1)"
    label variable is_core         "1 = core of division, 0 = overlap"

    sort geoid

    compress
    save  "$output/`stem'.dta", replace
    export delimited using "$output/`stem'.csv", replace

    display as result ///
        "report_10 (`lvl') done.  Exported " _N " units -> $output/`stem'.{csv,dta}"
}

display as result _n "report_10 done for: $levels."

cap log close
