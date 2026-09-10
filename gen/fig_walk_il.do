clear all
cap log close
local _logdate = subinstr("$S_DATE", " ", "", .)
log using "$log/fig_walk_il_${run_tag}_`_logdate'.log", replace text

* fig_walk_il.do -- draw the worked-example maps for Champaign County, IL and its ring of adjacent counties through the region_maps program.
* Called by: _master.do (generate phase, do_generate switch).
* Reads:  includes _inc/ensure_sp_dta.do and _inc/region_maps.do; via region_maps, $temp/tracts_sp.dta, $temp/counties_sp.dta and $temp/final_assignments_<run_tag>_tract.dta
* Writes: via region_maps, $figures/map_il_<hc|leiden|combined>_<run_tag>.png

include "$program/_inc/ensure_sp_dta.do"
include "$program/_inc/region_maps.do"

region_maps,                                                       ///
    focal("17019 17147 17041 17183")                               ///
    subtitle("Champaign County and adjacent counties, IL")         ///
    suffix("il")                                                   ///
    regionlabel("Champaign + ring")

display as result _n "fig_walk_il done."

cap log close
