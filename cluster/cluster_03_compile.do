clear all
cap log close
local _logdate = subinstr("$S_DATE", " ", "", .)
log using "$log/cluster_03_compile_${run_tag}_`_logdate'.log", replace text

* cluster_03_compile.do -- joins the HC and Leiden assignments at the chosen parameters into one file with a row per (unit, division) and the core/overlap flag.
* Called by: _master.do (cluster phase, $run_tag = baseline) and _inc/run_cluster.do (robust phase; $run_tag, $div_tag, $period_pflow_suffix and $hc_input set by the calling robust_* script).
* Reads:  $temp/chosen_params_<lvl>.dta; $temp/hc_assignments_<run_tag>_<lvl>.dta; $temp/leiden_assignments_<run_tag>_<lvl>.dta; $temp/division_lookup<tag>.dta; $temp/division_lookup_county<tag>.dta; includes _inc/ensure_county_div_lookup.do
* Writes: $temp/compiled_assignments_<run_tag>_<lvl>.dta
* Notes:  We abort if HC and Leiden do not cover the same (geoid, div_id) rows or if any row is absent from the division lookup, because a silent drop here would misalign the fusion step. A unit in an overlap zone still carries two rows, one per division.

include "$program/_inc/ensure_county_div_lookup.do"

foreach lvl of global levels {

    display as text _n "{hline 70}"
    display as text "===== cluster_03 LEVEL: `lvl' ====="
    display as text "{hline 70}"

    local _dtag  = cond("$div_tag" == "", "baseline", "$div_tag")
    local tagsuf = cond("`_dtag'" == "baseline", "", "_`_dtag'")
    if "`lvl'" == "tract" {
        local divlookup "$temp/division_lookup`tagsuf'.dta"
    }
    else {
        local divlookup "$temp/division_lookup_county`tagsuf'.dta"
    }

    use "$temp/chosen_params_`lvl'.dta", clear
    local k_nat      = hc_k_national[1]
    local res_chosen = leiden_resolution[1]
    display as text "Chosen hc_k_national = `k_nat'  leiden resolution = `res_chosen'"

    use "$temp/hc_assignments_${run_tag}_`lvl'.dta", clear
    rename hc_cluster hc_cluster_local
    keep geoid div_id hc_cluster_local

    tempfile hc
    save `hc'

    use "$temp/leiden_assignments_${run_tag}_`lvl'.dta", clear
    gen double _d = abs(resolution - `res_chosen')
    quietly summarize _d
    scalar _dmin = r(min)
    local res_dmin = r(min)
    if `res_dmin' > 0.5 * `res_chosen' {
        di as error "ERROR: nearest leiden resolution is far from chosen `res_chosen' (`lvl')."
        levelsof resolution, separate(", ")
        exit 459
    }
    keep if _d == _dmin
    scalar drop _dmin
    drop _d
    quietly count
    if r(N) == 0 {
        di as error "ERROR: no leiden_assignments rows near chosen resolution `res_chosen' (`lvl')."
        exit 459
    }
    rename leiden_cluster leiden_cluster_local
    keep geoid div_id leiden_cluster_local

    merge 1:1 geoid div_id using `hc'
    capture assert _merge == 3
    if _rc {
        quietly count if _merge != 3
        display as error "cluster_03: HC and Leiden assignments disagree on " ///
            r(N) " (geoid, div_id) rows (`lvl'); upstream node sets diverged."
        exit 459
    }
    drop _merge

    quietly count
    local _n_before = r(N)
    merge 1:1 geoid div_id using "`divlookup'", ///
        keepusing(is_core) keep(match) nogenerate
    quietly count
    if r(N) != `_n_before' {
        display as error "cluster_03: divlookup merge dropped " ///
            `=`_n_before' - r(N)' " assignment rows (`lvl'); assignments " ///
            "contain (geoid, div_id) pairs absent from the division lookup."
        exit 459
    }

    label variable hc_cluster_local     "HC cluster ID, local to division"
    label variable leiden_cluster_local "Leiden cluster ID, local to division"
    label variable is_core              "1 = core of window, 0 = overlap zone"

    order geoid div_id is_core hc_cluster_local leiden_cluster_local
    sort geoid div_id

    quietly count
    local nrows = r(N)
    quietly count if is_core == 0
    local n_overlap = r(N)

    bysort geoid: gen long _ndivs = _N
    quietly count if _ndivs > 1
    local n_dup_rows = r(N)
    drop _ndivs
    sort geoid div_id

    display as text _n "Compiled assignments (`lvl'):"
    display as text "  Total rows:                       " %8.0fc `nrows'
    display as text "  Rows in overlap zones:            " %8.0fc `n_overlap'
    display as text "  Units assigned to >1 division:    " %8.0fc `n_dup_rows'

    compress
    save "$temp/compiled_assignments_${run_tag}_`lvl'.dta", replace
    display as result "cluster_03 (`lvl') done."
}

display as result _n "cluster_03 done for $levels."

cap log close
