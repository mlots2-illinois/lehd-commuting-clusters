clear all
cap log close
local _logdate = subinstr("$S_DATE", " ", "", .)
log using "$log/fig_flow_map_${run_tag}_`_logdate'.log", replace text

* fig_flow_map.do -- map the dominant outbound commuting flow of every tract in Champaign-Urbana, IL (Champaign County, FIPS 17019) as arrows over the tract outlines.
* Called by: _master.do (generate phase, do_generate switch).
* Reads:  includes _inc/ensure_sp_dta.do; $temp/tracts_sp.dta; $temp/tracts_sp_shp.dta; $temp/centroids_tract.dta; $temp/pflow_tract_<period_pflow_suffix>.dta
* Writes: $figures/flow_map_il_<run_tag>.png
* Notes:  We keep only the single largest destination per origin tract and bin those flows into terciles of P_ij, so the arrow weight reads as flow strength. The bounding box trims the county to the urban core, and one outlying tract (17019010604) is dropped by geoid.

include "$program/_inc/ensure_sp_dta.do"
ensure_sp_dta, level(tract)
local _ok_tract = r(ok)
ensure_sp_dta, level(county)
local _ok_county = r(ok)
if !`_ok_tract' | !`_ok_county' {
    di as text "fig_flow_map: shapefile(s) unavailable; skipping."
    cap log close
    exit 0
}

local focal "17019"
local pflow "$temp/pflow_tract_${period_pflow_suffix}.dta"

capture confirm file "`pflow'"
if _rc {
    di as text "fig_flow_map: `pflow' not found (run build_01 at tract level); skipping."
    cap log close
    exit 0
}
capture confirm file "$temp/centroids_tract.dta"
if _rc {
    di as text "fig_flow_map: $temp/centroids_tract.dta not found (run build_02 at tract level); skipping."
    cap log close
    exit 0
}

use _ID GEOID _CX _CY using "$temp/tracts_sp.dta", clear
rename GEOID geoid
merge 1:1 geoid using "$temp/centroids_tract.dta", ///
    keep(master match) keepusing(lat lon pop) nogenerate

local cu_latlo  40.05
local cu_lathi  40.17
local cu_lonlo -88.31
local cu_lonhi -88.13
local drop_tracts "17019010604"
keep if substr(geoid, 1, 5) == "`focal'" ///
    & inrange(lat, `cu_latlo', `cu_lathi') & inrange(lon, `cu_lonlo', `cu_lonhi') ///
    & !inlist(geoid, "17019010604")
replace pop = 0 if mi(pop)
quietly count
display as text "Focal (Champaign-Urbana proper) tracts: " r(N)

quietly summarize _ID
local tr_idmin = r(min)
local tr_idmax = r(max)
preserve
keep _ID
tempfile focal_tr_ids
save `focal_tr_ids'
restore
drop _ID
tempfile cent
save `cent'

use geo_i geo_j P_ij using "`pflow'", clear
rename geo_i geoid
merge m:1 geoid using `cent', keep(match) keepusing(_CX _CY) nogenerate
rename (geoid _CX _CY) (geo_i cx_i cy_i)
rename geo_j geoid
merge m:1 geoid using `cent', keep(match) keepusing(_CX _CY) nogenerate
rename (geoid _CX _CY) (geo_j cx_j cy_j)

quietly count
if r(N) == 0 {
    di as text "fig_flow_map: no within-region flow pairs found; check the focal county list. Skipping."
    cap log close
    exit 0
}
display as text "Within-region flow pairs: " r(N)

drop if geo_i == geo_j
bysort geo_i (P_ij): keep if _n == _N
quietly count
display as text "Dominant-flow lines (one per origin tract): " r(N)

xtile pbin = P_ij, nq(3)
tempfile flows
save `flows'

use `cent', clear
rename (_CX _CY pop) (ncx ncy npop)
keep ncx ncy npop
append using `flows'
tempfile layers
save `layers'

use _ID _X _Y using "$temp/tracts_sp_shp.dta" ///
    if inrange(_ID, `tr_idmin', `tr_idmax'), clear
merge m:1 _ID using `focal_tr_ids', keep(match) nogenerate
keep _ID _X _Y
rename (_X _Y) (tcx tcy)
tempfile tract_outline
save `tract_outline'

use `tract_outline', clear
append using `layers'

twoway ///
    (line tcy tcx, cmissing(n) lcolor(gs11) lwidth(vthin)) ///
    (pcarrow cy_i cx_i cy_j cx_j if pbin == 1, ///
    lcolor("32 64 128 %35") lwidth(vthin) mcolor("32 64 128 %35") msize(1.2)) ///
    (pcarrow cy_i cx_i cy_j cx_j if pbin == 2, ///
    lcolor("32 64 128 %60") lwidth(thin) mcolor("32 64 128 %60") msize(1.4)) ///
    (pcarrow cy_i cx_i cy_j cx_j if pbin == 3, ///
    lcolor("32 64 128 %90") lwidth(medthin) mcolor("32 64 128 %90") msize(1.6)) ///
    (scatter ncy ncx [aw=npop], ///
    msymbol(circle) mcolor("139 0 0 %55") mlcolor(white) mlwidth(vvthin)) ///
    , aspectratio(1) legend(off) ///
    xtitle("") ytitle("") xlabel(none) ylabel(none) ///
    xscale(off) yscale(off) ///
    graphregion(color(white)) plotregion(lstyle(none))

graph export "$figures/flow_map_il_${run_tag}.png", replace width(2400)
display as result "  -> flow_map_il_${run_tag}.png"

display as result _n "fig_flow_map done."

cap log close
