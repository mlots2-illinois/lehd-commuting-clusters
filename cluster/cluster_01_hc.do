clear all
cap log close
local _logdate = subinstr("$S_DATE", " ", "", .)
log using "$log/cluster_01_hc_${run_tag}_`_logdate'.log", replace text

* cluster_01_hc.do -- builds one dissimilarity matrix per division and cuts a hierarchical clustering (HC) of each at the k allocated from the national target.
* Called by: _master.do (cluster phase, $run_tag = baseline) and _inc/run_cluster.do (robust phase; $run_tag, $div_tag, $period_pflow_suffix and $hc_input set by the calling robust_* script).
* Reads:  $temp/chosen_params_<lvl>.dta; $temp/pflow_<lvl>_<period_pflow_suffix>.dta; $temp/division_lookup<tag>.dta; $temp/division_lookup_county<tag>.dta; includes _inc/ensure_county_div_lookup.do, _inc/build_div_dissim.do
* Writes: via _py/build_dissim.py: $hc_input/<lvl>_<run_tag>/dissim_div<div>.csv, $temp/_c01_counts_<run_tag>_<lvl>.csv; via _py/hc_run.py: $temp/hc_output_<run_tag>_<lvl>.csv; $temp/hc_assignments_<run_tag>_<lvl>.dta
* Notes:  The per-division k is round(n_units / (N / k_national)), floored at 2 and capped at n_units, the rule calib_03 also uses. We set every BLAS thread pool to one before the Python steps so that the dissimilarity build and the HC run do not oversubscribe the cores. The node-set check runs for the baseline only, since the robust reruns share its geometry.

include "$program/_inc/ensure_county_div_lookup.do"
include "$program/_inc/build_div_dissim.do"

local pyscript "$program/_py/hc_run.py"

capture mkdir "$hc_input"

foreach lvl of global levels {

    display as text _n "{hline 70}"
    display as text "===== cluster_01 LEVEL: `lvl' ====="
    display as text "{hline 70}"

    local _dtag  = cond("$div_tag" == "", "baseline", "$div_tag")
    local tagsuf = cond("`_dtag'" == "baseline", "", "_`_dtag'")
    if "`lvl'" == "tract" {
        local divlookup "$temp/division_lookup`tagsuf'.dta"
    }
    else {
        local divlookup "$temp/division_lookup_county`tagsuf'.dta"
    }
    local pflow_file "$temp/pflow_`lvl'_${period_pflow_suffix}.dta"
    local params     "$temp/chosen_params_`lvl'.dta"

    local hc_dir "$hc_input/`lvl'_${run_tag}"
    capture mkdir "`hc_dir'"

    local out_csv "$temp/hc_output_${run_tag}_`lvl'.csv"

    use "`params'", clear
    local alpha     = alpha[1]
    local cap_km    = cap_km[1]
    local k_nat     = hc_k_national[1]

    display as text "alpha = `alpha'  cap_km = `cap_km'  hc_k_national = `k_nat'"

    use geo_i geo_j P_ij using "`pflow_file'", clear
    tempfile pflow_all
    save `pflow_all'

    use div_id geoid lat lon using "`divlookup'", clear

    quietly {
        bysort geoid (div_id): gen byte _ufirst = (_n == 1)
        count if _ufirst
        local total_n = r(N)
        drop _ufirst
    }
    local avg_size = `total_n' / `k_nat'
    display as text "Total `lvl' units: `total_n'   target avg cluster size: " %6.2f `avg_size'

    quietly levelsof div_id, local(divs)
    local ndivs : word count `divs'
    display as text "Writing dissim CSVs for `ndivs' divisions ..."

    if "${run_tag}" == "baseline" {
        preserve
        keep geoid
        duplicates drop
        gen byte _in_geom = 1
        tempfile _geomnodes
        save `_geomnodes'

        use geo_i using `pflow_all', clear
        rename geo_i geoid
        tempfile _fi
        save `_fi'
        use geo_j using `pflow_all', clear
        rename geo_j geoid
        append using `_fi'
        keep geoid
        duplicates drop
        gen byte _in_flow = 1

        merge 1:1 geoid using `_geomnodes', nogenerate
        quietly count if _in_flow == 1 & mi(_in_geom)
        local _flow_only = r(N)
        quietly count if _in_geom == 1 & mi(_in_flow)
        local _geom_only = r(N)
        display as text "  [node-set] flow-only (in graph, no centroid/division): `_flow_only'"
        display as text "  [node-set] geometry-only (centroid, no flow edges):    `_geom_only'"
        restore
    }

    local pyenvpfx "OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1"
    local counts_csv "$temp/_c01_counts_${run_tag}_`lvl'.csv"
    capture erase "`counts_csv'"
    shell `pyenvpfx' "$py" "$program/_py/build_dissim.py" ///
        "`divlookup'" "`pflow_file'" "`hc_dir'" `alpha' `cap_km' "$divfmt" "`counts_csv'"
    capture confirm file "`counts_csv'"
    if _rc {
        display as error "cluster_01: build_dissim.py produced no counts (`counts_csv'); check $py env and the log."
        exit 601
    }
    capture frame drop _c01_counts
    frame create _c01_counts
    frame _c01_counts: import delimited using "`counts_csv'", varnames(1) clear

    local k_map_str ""
    frame _c01_counts {
        local _nc = _N
        forvalues i = 1/`_nc' {
            local d     = div_id[`i']
            local N     = n_units[`i']
            local k_div = round(`N' / `avg_size')
            if `k_div' < 2   local k_div 2
            if `k_div' > `N' local k_div `N'
            if "`k_map_str'" == "" {
                local k_map_str "`d':`k_div'"
            }
            else {
                local k_map_str "`k_map_str',`d':`k_div'"
            }
        }
    }
    frame drop _c01_counts

    display as text _n "Calling: $py `pyscript' `hc_dir' --k-map ..."
    capture erase "`out_csv'"
    shell `pyenvpfx' "$py" "`pyscript'" "`hc_dir'" "`out_csv'" --k-map "`k_map_str'"
    capture confirm file "`out_csv'"
    if _rc {
        display as error "cluster_01: Python step produced no output (`out_csv'); check $py env and the log."
        exit 601
    }

    import delimited using "`out_csv'", varnames(1) stringcols(2) clear
    confirm string variable geoid
    assert strlen(geoid) == cond("`lvl'" == "tract", 11, 5)
    rename cut_value k_requested

    label variable div_id      "Division ID"
    label variable geoid       "Geographic FIPS (`lvl')"
    label variable k_requested "Per-division k requested (target maxclust)"
    label variable hc_cluster  "Cluster ID within division"

    order geoid div_id k_requested hc_cluster
    sort geoid div_id
    compress
    save "$temp/hc_assignments_${run_tag}_`lvl'.dta", replace

    display as result _n ///
        "cluster_01 (`lvl') done. Saved hc_assignments_${run_tag}_`lvl'.dta: " _N " rows."
}

display as result _n "cluster_01 done for $levels."

cap log close
