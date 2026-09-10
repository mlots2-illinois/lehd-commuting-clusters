clear all
cap log close
local _logdate = subinstr("$S_DATE", " ", "", .)
log using "$log/fig_walk_chicago_${run_tag}_`_logdate'.log", replace text

* fig_walk_chicago.do -- draw the worked-example maps for the six-county Chicago metropolitan core through the region_maps program.
* Called by: _master.do (generate phase, do_generate switch).
* Reads:  includes _inc/ensure_sp_dta.do and _inc/region_maps.do; via region_maps, $temp/tracts_sp.dta, $temp/counties_sp.dta and $temp/final_assignments_<run_tag>_tract.dta
* Writes: via region_maps, $figures/map_chicago_<hc|leiden|combined>_<run_tag>.png
* Notes:  Cook County alone holds about 1,300 tracts, so this example shows how fine the clusters become in a dense urban core.

include "$program/_inc/ensure_sp_dta.do"
include "$program/_inc/region_maps.do"

region_maps,                                                       ///
    focal("17031 17043 17089 17097 17111 17197")                   ///
    subtitle("Chicago metropolitan core (six-county Chicagoland), IL") ///
    suffix("chicago")                                              ///
    regionlabel("Chicago metro")                                   ///
    note2("Cook County (Chicago) holds ~1,300 tracts; clusters are finest in the dense urban core.")

display as result _n "fig_walk_chicago done."

cap log close
