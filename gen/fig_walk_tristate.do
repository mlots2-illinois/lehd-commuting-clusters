clear all
cap log close
local _logdate = subinstr("$S_DATE", " ", "", .)
log using "$log/fig_walk_tristate_${run_tag}_`_logdate'.log", replace text

* fig_walk_tristate.do -- draw the worked-example maps for the Quad Cities region of eastern Iowa and western Illinois through the region_maps program.
* Called by: _master.do (generate phase, do_generate switch).
* Reads:  includes _inc/ensure_sp_dta.do and _inc/region_maps.do; via region_maps, $temp/tracts_sp.dta, $temp/counties_sp.dta and $temp/final_assignments_<run_tag>_tract.dta
* Writes: via region_maps, $figures/map_tristate_<hc|leiden|combined>_<run_tag>.png
* Notes:  This example shows a commuting cluster that spans the Iowa-Illinois line at the Mississippi River.

include "$program/_inc/ensure_sp_dta.do"
include "$program/_inc/region_maps.do"

region_maps,                                                                  ///
    focal("19031 19045 19163 19139 19115 17195 17161 17073 17131")            ///
    subtitle("Quad Cities region: Eastern Iowa and Western Illinois")         ///
    suffix("tristate")                                                        ///
    regionlabel("Quad Cities region")                                         ///
    note2("The Iowa-Illinois state line follows the Mississippi River between Scott (IA) and Rock Island (IL).")

display as result _n "fig_walk_tristate done."

cap log close
