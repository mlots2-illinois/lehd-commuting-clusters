clear all
cap log close
local _logdate = subinstr("$S_DATE", " ", "", .)
log using "$log/fig_map_regions_${run_tag}_`_logdate'.log", replace text

* fig_map_regions.do -- draw the tract partition for each of the four Census regions through the conus_map program.
* Called by: _master.do (generate phase, do_generate switch).
* Reads:  includes _inc/ensure_sp_dta.do, _inc/conus_map.do and _inc/census_regions.do; via conus_map, $temp/tracts_sp.dta and $temp/final_assignments_<run_tag>_tract.dta
* Writes: via conus_map, $figures/map_<region>_hc_<run_tag>.pdf and $figures/map_<region>_leiden_<run_tag>.pdf
* Notes:  The region-to-state lists come from $states_<region> in _inc/census_regions.do.

include "$program/_inc/ensure_sp_dta.do"
include "$program/_inc/conus_map.do"
include "$program/_inc/census_regions.do"

foreach r of global region_keys {
    display as text _n "{hline 70}"
    display as text "fig_map_regions: `r' (states ${states_`r'})"
    display as text "{hline 70}"
    conus_map, suffix("`r'") keepstates("${states_`r'}") level(tract)
}

display as result _n "fig_map_regions done.  Regional maps written to $figures/."

cap log close
