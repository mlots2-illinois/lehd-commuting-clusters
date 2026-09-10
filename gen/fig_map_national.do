clear all
cap log close
local _logdate = subinstr("$S_DATE", " ", "", .)
log using "$log/fig_map_national_${run_tag}_`_logdate'.log", replace text

* fig_map_national.do -- draw the conterminous US map of the HC and Leiden partitions through the conus_map program.
* Called by: _master.do (generate phase, do_generate switch).
* Reads:  includes _inc/ensure_sp_dta.do and _inc/conus_map.do; via conus_map, $temp/counties_sp.dta and $temp/final_assignments_<run_tag>_county.dta
* Writes: via conus_map, $figures/map_national_hc_<run_tag>.pdf and $figures/map_national_leiden_<run_tag>.pdf
* Notes:  No level() is passed, so conus_map draws its default, the county partition; the tract partition is mapped by region in fig_map_regions.do.

include "$program/_inc/ensure_sp_dta.do"
include "$program/_inc/conus_map.do"

conus_map, suffix("national")

display as result _n "fig_map_national done.  National maps written to $figures/."

cap log close
