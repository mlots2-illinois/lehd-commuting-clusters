clear all
cap log close
local _logdate = subinstr("$S_DATE", " ", "", .)
log using "$log/cluster_06_harmonize_${run_tag}_`_logdate'.log", replace text

* cluster_06_harmonize.do -- replaces the division-local HC and Leiden cluster IDs in the reconciled file with the global IDs from the fusion maps.
* Called by: _master.do (cluster phase, $run_tag = baseline) and _inc/run_cluster.do (robust phase; $run_tag, $div_tag, $period_pflow_suffix and $hc_input set by the calling robust_* script).
* Reads:  $temp/reconciled_assignments_<run_tag>_<lvl>.dta; $temp/fusion_map_hc_<run_tag>_<lvl>.dta; $temp/fusion_map_leiden_<run_tag>_<lvl>.dta
* Writes: $temp/harmonized_assignments_<run_tag>_<lvl>.dta
* Notes:  A unit whose (division, local cluster) pair is missing from a fusion map receives a fresh ID above the current maximum. That path is a warning, not an expected case, since cluster_04 maps every local cluster.

foreach lvl of global levels {

    display as text _n "{hline 70}"
    display as text "===== cluster_06 LEVEL: `lvl' ====="
    display as text "{hline 70}"

    use "$temp/reconciled_assignments_${run_tag}_`lvl'.dta", clear

    rename hc_cluster_local local_cluster
    merge m:1 div_id local_cluster using "$temp/fusion_map_hc_${run_tag}_`lvl'.dta", ///
        keep(master match) keepusing(global_cluster) nogenerate
    rename global_cluster hc_cluster
    rename local_cluster  hc_cluster_local

    rename leiden_cluster_local local_cluster
    merge m:1 div_id local_cluster using "$temp/fusion_map_leiden_${run_tag}_`lvl'.dta", ///
        keep(master match) keepusing(global_cluster) nogenerate
    rename global_cluster leiden_cluster
    rename local_cluster  leiden_cluster_local

    foreach algo in hc leiden {
        quietly count if mi(`algo'_cluster)
        if r(N) > 0 {
            di as error "WARNING: `r(N)' units missing `algo' global ID; assigning fresh unique fallback IDs"
            quietly summarize `algo'_cluster
            local _base = cond(r(N) > 0 & r(max) < ., r(max), 0)
            egen long _fb = group(div_id `algo'_cluster_local) if mi(`algo'_cluster)
            quietly replace `algo'_cluster = `_base' + _fb if mi(`algo'_cluster)
            drop _fb
        }
    }

    drop hc_cluster_local leiden_cluster_local

    label variable hc_cluster     "Globally unique HC cluster ID (post-fusion)"
    label variable leiden_cluster "Globally unique Leiden cluster ID (post-fusion)"

    order geoid div_id is_core hc_cluster leiden_cluster
    sort geoid

    egen byte _tag_hc = tag(hc_cluster)
    quietly count if _tag_hc
    local hc_n = r(N)
    egen byte _tag_le = tag(leiden_cluster)
    quietly count if _tag_le
    local le_n = r(N)
    drop _tag_hc _tag_le
    display as text _n "Harmonized cluster counts (`lvl'):"
    display as text "  HC clusters:     " %6.0f `hc_n'
    display as text "  Leiden clusters: " %6.0f `le_n'
    display as text "  Units:           " %6.0f _N

    compress
    save "$temp/harmonized_assignments_${run_tag}_`lvl'.dta", replace
    display as result "cluster_06 (`lvl') done."
}

display as result _n "cluster_06 done for $levels."

cap log close
