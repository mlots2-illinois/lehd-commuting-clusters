clear all
cap log close
local _logdate = subinstr("$S_DATE", " ", "", .)
log using "$log/fig_walk_omaha_${run_tag}_`_logdate'.log", replace text

* fig_walk_omaha.do -- draw the worked-example maps for Omaha-Council Bluffs, spanning eastern Nebraska and western Iowa, through the region_maps program.
* Called by: _master.do (generate phase, do_generate switch).
* Reads:  includes _inc/ensure_sp_dta.do and _inc/region_maps.do; via region_maps, $temp/tracts_sp.dta, $temp/counties_sp.dta and $temp/final_assignments_<run_tag>_tract.dta
* Writes: via region_maps, $figures/map_omaha_<hc|leiden|combined>_<run_tag>.png
* Notes:  This example shows a commuting cluster that spans the Nebraska-Iowa line at the Missouri River.

include "$program/_inc/ensure_sp_dta.do"
include "$program/_inc/region_maps.do"

region_maps,                                                                  ///
    focal("31055 31153 31177 31155 31025 19155 19129 19085")                  ///
    subtitle("Eastern Nebraska + Western Iowa") ///
    suffix("omaha")                                                           ///
    regionlabel("Omaha–Council Bluffs")                                       ///
    note2("The Nebraska–Iowa state line follows the Missouri River between Douglas (NE) and Pottawattamie (IA).")

display as result _n "fig_walk_omaha done."

cap log close
