clear all
cap log close
local _logdate = subinstr("$S_DATE", " ", "", .)
log using "$log/report_04_validation_${run_tag}_`_logdate'.log", replace text

* report_04_validation.do -- run four data-integrity checks on the final assignments (coverage, no missing cluster IDs, no duplicate units, plausible cluster counts) and write a PASS/FAIL report per geography.
* Called by: _master.do (generate phase, do_generate switch).
* Reads:  $temp/centroids_<lvl>.dta; $temp/final_assignments_<run_tag>_<lvl>.dta
* Writes: $output/validation_report_<run_tag>_<lvl>.txt
* Notes:  Any failed check stops the run with error 459 before the export phase, so a deliverable is never written from a partition that failed integrity.

local total_fails 0

foreach lvl of global levels {

    display as text _n "{hline 70}"
    display as text "===== report_04 LEVEL: `lvl' ====="
    display as text "{hline 70}"

    capture confirm file "$temp/centroids_`lvl'.dta"
    if _rc {
        display as text "report_04: $temp/centroids_`lvl'.dta not found; skipping level `lvl'."
        continue
    }
    capture confirm file "$temp/final_assignments_${run_tag}_`lvl'.dta"
    if _rc {
        display as text "report_04: no final_assignments for `lvl'; skipping level."
        continue
    }

    local fail_count 0
    local report "$output/validation_report_${run_tag}_`lvl'.txt"
    capture file close vlog
    file open vlog using "`report'", write replace

    file write vlog "Validation report for run_tag = ${run_tag}, level = `lvl'" _n
    file write vlog "Run date: " "`c(current_date)' `c(current_time)'" _n _n

    use geoid using "$temp/centroids_`lvl'.dta", clear
    quietly count
    local n_universe = r(N)

    merge 1:1 geoid using "$temp/final_assignments_${run_tag}_`lvl'.dta", ///
        keepusing(hc_cluster) generate(_m)

    quietly count if _m == 1
    local n_only_centroids = r(N)
    quietly count if _m == 2
    local n_only_assigned  = r(N)
    quietly count if _m == 3
    local n_matched        = r(N)

    local check1 = (`n_only_centroids' == 0) & (`n_only_assigned' == 0)
    local status = cond(`check1', "PASS", "FAIL")
    file write vlog "[`status'] 1. Coverage: " ///
        "centroids=`n_universe'  matched=`n_matched'  " ///
        "centroids-only=`n_only_centroids'  assigned-only=`n_only_assigned'" _n
    display as text "[`status'] coverage: " `n_matched' " of " `n_universe'
    if !`check1' local ++fail_count

    use "$temp/final_assignments_${run_tag}_`lvl'.dta", clear
    quietly count if mi(hc_cluster)
    local n_mi_hc = r(N)
    quietly count if mi(leiden_cluster)
    local n_mi_le = r(N)

    local check2 = (`n_mi_hc' == 0) & (`n_mi_le' == 0)
    local status = cond(`check2', "PASS", "FAIL")
    file write vlog "[`status'] 2. Cluster ID coverage: " ///
        "missing HC=`n_mi_hc'  missing Leiden=`n_mi_le'" _n
    display as text "[`status'] no missing cluster IDs"
    if !`check2' local ++fail_count

    quietly duplicates report geoid
    local n_dup = r(N) - r(unique_value)

    local check3 = (`n_dup' == 0)
    local status = cond(`check3', "PASS", "FAIL")
    file write vlog "[`status'] 3. Duplicate geoid rows: `n_dup'" _n
    display as text "[`status'] no duplicate unit rows"
    if !`check3' local ++fail_count

    egen byte _tag_hc = tag(hc_cluster)
    quietly count if _tag_hc
    local n_hc_clusters = r(N)
    egen byte _tag_le = tag(leiden_cluster)
    quietly count if _tag_le
    local n_le_clusters = r(N)
    drop _tag_hc _tag_le
    quietly count
    local n_units = r(N)

    local check4 = (`n_hc_clusters' > 0)         ///
        & (`n_hc_clusters' < `n_units') ///
        & (`n_le_clusters' > 0)         ///
        & (`n_le_clusters' < `n_units')
    local status = cond(`check4', "PASS", "FAIL")
    file write vlog "[`status'] 4. Cluster counts: " ///
        "HC=`n_hc_clusters'  Leiden=`n_le_clusters'  units=`n_units'" _n
    display as text "[`status'] cluster count plausibility"
    if !`check4' local ++fail_count

    file write vlog _n
    if `fail_count' == 0 {
        file write vlog "OVERALL (`lvl'): PASS" _n
        display as result "Validation (`lvl'): PASS"
    }
    else {
        file write vlog "OVERALL (`lvl'): FAIL (`fail_count' check(s) failed)" _n
        display as error "Validation (`lvl') FAILED: `fail_count' check(s) failed"
    }
    file close vlog
    display as text "Wrote: `report'"

    local total_fails = `total_fails' + `fail_count'
}

display as result _n "report_04 done (data-integrity checks only; CBSA preservation removed)."
cap log close

if `total_fails' > 0 {
    display as error "report_04: `total_fails' validation check(s) FAILED; see the validation_report_*.txt files."
    exit 459
}
