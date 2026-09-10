clear all
cap log close
local _logdate = subinstr("$S_DATE", " ", "", .)
log using "$log/calib_05_validity_${run_tag}_`_logdate'.log", replace text

* calib_05_validity.do -- adds bootstrap stability and null-model checks for the HC cut, the validity evidence that the monotone selectors in calib_03 cannot supply.
* Called by: _master.do (calib phase, when $do_calib_validity; after calib_04; expects $n_sample_divs, $k_grid_<lvl>, $hs_nboot, $nt_nperm).
* Reads:  $temp/division_lookup.dta; $temp/division_lookup_county.dta; $temp/pflow_<lvl>_<period_pflow_suffix>.dta; includes _inc/calib_params.do, _inc/sample_divisions.do
* Writes: via _py/hc_stability.py: $tables/hc_stability_<lvl>.csv; via _py/null_structure_test.py: $tables/null_structure_<lvl>.csv; via _py/calib_hc_k_view.py: $figures/calib_hc_k_view_<lvl>.png (refreshed)
* Notes:  Self-containment and the within-cluster spread guard are both monotone in k, so neither can certify the chosen value. hc_stability.py runs the clusterwise bootstrap of Hennig (2007), the HC analogue of the cross-seed reproducibility calib_02 reports for Leiden. null_structure_test.py compares observed self-containment with a geography-only partition and a permutation null (Hennig and Lin 2015). Silhouette and modularity over the same k grid come from calib_03 with PQ_FULL=1.

include "$program/_inc/calib_params.do"
include "$program/_inc/sample_divisions.do"

* Defaults let the script run on its own; _master.do sets B = 50 and R = 99 explicitly.
local nboot  = cond("$hs_nboot"  == "", 50, $hs_nboot)
local nperm  = cond("$nt_nperm"  == "", 99, $nt_nperm)

foreach lvl of global levels {

    display as text _n "{hline 70}"
    display as text "===== calib_05 LEVEL: `lvl' ====="
    display as text "{hline 70}"

    if "`lvl'" == "tract" {
        local divlookup "$temp/division_lookup.dta"
        local k_grid    "$k_grid_tract"
    }
    else {
        local divlookup "$temp/division_lookup_county.dta"
        local k_grid    "$k_grid_county"
    }
    local alpha     = ``lvl'_alpha'
    local cap_km    = ``lvl'_cap_km'
    local chosen_k  = ``lvl'_hc_k_national'
    local pflow     "$temp/pflow_`lvl'_${period_pflow_suffix}.dta"

    local _ok 1
    foreach f in "`divlookup'" "`pflow'" {
        capture confirm file "`f'"
        if _rc {
            display as error "  missing `f'; skipping `lvl'."
            local _ok 0
        }
    }
    if !`_ok' continue

    * One row per unit, core division preferred, so the sampled divisions match those of calib_03.
    use div_id geoid is_core using "`divlookup'", clear
    gen byte _neg_core = -is_core
    bysort geoid (_neg_core div_id): keep if _n == 1
    sample_divisions, n($n_sample_divs)
    local samp_csv = subinstr(trim(itrim("`r(sample_divs)'")), " ", ",", .)

    local kgrid_csv = subinstr(trim(itrim("`k_grid'")), " ", ",", .)

    * ---- clusterwise bootstrap stability -------------------------------
    local st_csv "$tables/hc_stability_`lvl'.csv"
    capture erase "`st_csv'"
    display as text "  bootstrap stability (B=`nboot') over k grid ..."
    shell HS_KGRID="`kgrid_csv'" HS_SAMPLE_DIVS="`samp_csv'" HS_NBOOT=`nboot' ///
        "$py" "$program/_py/hc_stability.py" "`lvl'" "`divlookup'" "`pflow'" ///
        "`cap_km'" "`alpha'" "$seed" "`st_csv'"
    capture confirm file "`st_csv'"
    if _rc  display as error "  hc_stability.py produced no output for `lvl'."
    else    display as result "  -> `st_csv'"

    * ---- geography-only and permutation null ---------------------------
    local nt_csv "$tables/null_structure_`lvl'.csv"
    capture erase "`nt_csv'"
    display as text "  permutation null (R=`nperm') over k grid ..."
    shell NT_KGRID="`kgrid_csv'" NT_SAMPLE_DIVS="`samp_csv'" NT_NPERM=`nperm' ///
        "$py" "$program/_py/null_structure_test.py" "`lvl'" "`divlookup'" ///
        "`pflow'" "`cap_km'" "`alpha'" "$seed" "`nt_csv'"
    capture confirm file "`nt_csv'"
    if _rc  display as error "  null_structure_test.py produced no output for `lvl'."
    else    display as result "  -> `nt_csv'"

    * ---- refresh the k view with the validity and stability panels -----
    * calib_hc_k_view.py reads the stability and null CSVs when they exist, so the k figure gains its two extra panels only after this script.
    capture confirm file "$program/_py/calib_hc_k_view.py"
    if !_rc {
        shell "$py" "$program/_py/calib_hc_k_view.py" "$tables" "$figures" ///
            "`lvl'" `chosen_k'
    }
}

display as result _n "calib_05 done. Stability and null tables in $tables/."

cap log close
