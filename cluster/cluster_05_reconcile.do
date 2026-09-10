clear all
cap log close
local _logdate = subinstr("$S_DATE", " ", "", .)
log using "$log/cluster_05_reconcile_${run_tag}_`_logdate'.log", replace text

* cluster_05_reconcile.do -- resolves each overlap unit to the one division whose members it commutes with most, leaving one row per unit.
* Called by: _master.do (cluster phase, $run_tag = baseline) and _inc/run_cluster.do (robust phase; $run_tag, $div_tag, $period_pflow_suffix and $hc_input set by the calling robust_* script).
* Reads:  $temp/compiled_assignments_<run_tag>_<lvl>.dta; $temp/pflow_<lvl>_<period_pflow_suffix>.dta; $temp/division_lookup<tag>.dta; $temp/division_lookup_county<tag>.dta; includes _inc/ensure_county_div_lookup.do
* Writes: $temp/reconciled_assignments_<run_tag>_<lvl>.dta
* Notes:  A candidate division's score is the sum of P_ij between the unit and every unit in that division, so we assign by observed commuting rather than by the is_core flag from build_03. Ties, and units with no positive flow to any candidate, fall to the lowest div_id.

include "$program/_inc/ensure_county_div_lookup.do"

foreach lvl of global levels {

    display as text _n "{hline 70}"
    display as text "===== cluster_05 LEVEL: `lvl' ====="
    display as text "{hline 70}"

    local _dtag  = cond("$div_tag" == "", "baseline", "$div_tag")
    local tagsuf = cond("`_dtag'" == "baseline", "", "_`_dtag'")
    if "`lvl'" == "tract" {
        local divlookup "$temp/division_lookup`tagsuf'.dta"
    }
    else {
        local divlookup "$temp/division_lookup_county`tagsuf'.dta"
    }
    local pflow    "$temp/pflow_`lvl'_${period_pflow_suffix}.dta"
    local src      "$temp/compiled_assignments_${run_tag}_`lvl'.dta"
    local outfile  "$temp/reconciled_assignments_${run_tag}_`lvl'.dta"

    use "`src'", clear
    bysort geoid: gen long _ndivs = _N
    quietly count if _ndivs > 1
    local has_overlap = (r(N) > 0)

    if !`has_overlap' {
        display as text "No overlap units (`lvl'); passing compiled through."
        drop _ndivs
        isid geoid
        compress
        save "`outfile'", replace
        display as result "cluster_05 (`lvl') done. " _N " unique units."
        continue
    }

    keep if _ndivs > 1
    keep geoid div_id
    rename geoid t_geoid
    rename div_id candidate_div

    tempfile candidates
    save `candidates'

    quietly count
    display as text "Overlap candidate rows: " _N

    keep t_geoid
    duplicates drop
    tempfile overlap_units
    save `overlap_units'

    * pflow stores each pair once with geo_i < geo_j, so we read it twice to see the unit on either side.
    use geo_i geo_j P_ij using "`pflow'", clear
    rename geo_i t_geoid
    rename geo_j partner
    merge m:1 t_geoid using `overlap_units', keep(match) nogenerate
    tempfile dirA
    save `dirA'

    use geo_i geo_j P_ij using "`pflow'", clear
    rename geo_j t_geoid
    rename geo_i partner
    merge m:1 t_geoid using `overlap_units', keep(match) nogenerate
    append using `dirA'

    preserve
    use geoid div_id using "`divlookup'", clear
    rename geoid partner
    rename div_id partner_div
    tempfile pdivs
    save `pdivs'
    restore

    * A partner that sits in an overlap zone belongs to two divisions and contributes its flow to both candidate scores.
    joinby partner using `pdivs', unmatched(none)

    joinby t_geoid using `candidates', unmatched(none)
    keep if partner_div == candidate_div

    gcollapse (sum) score = P_ij, by(t_geoid candidate_div)

    * Highest score wins; on a tie the lower div_id, so the result does not depend on sort order.
    gsort t_geoid -score candidate_div
    by t_geoid: keep if _n == 1
    rename t_geoid geoid
    rename candidate_div assigned_div
    keep geoid assigned_div

    tempfile picked
    save `picked'
    display as text "Overlap units with a positive flow score: " _N

    use "`src'", clear
    merge m:1 geoid using `picked', keep(master match) nogenerate

    bysort geoid: gen long _ndivs = _N
    gen byte _drop_winner   = (_ndivs > 1) & !mi(assigned_div) & (div_id != assigned_div)
    * Units with no positive flow to any candidate keep their lowest div_id row.
    bysort geoid (div_id): gen byte _drop_fallback = (_ndivs > 1) & mi(assigned_div) & (_n > 1)

    drop if _drop_winner | _drop_fallback
    drop _drop_winner _drop_fallback assigned_div _ndivs

    isid geoid

    label data "Reconciled cluster assignments (one row per unit)"
    compress
    save "`outfile'", replace
    display as result "cluster_05 (`lvl') done. " _N " unique units."
}

display as result _n "cluster_05 done for $levels."

cap log close
