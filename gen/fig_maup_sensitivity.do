clear all
cap log close

* fig_maup_sensitivity.do -- plot the earnings-on-BA slope with its 95 percent interval at each areal unit on the assess_09_maup.do ladder, coarsest to finest.
* Called by: _master.do (generate phase, do_generate switch).
* Reads:  $tables/maup_<run_tag>.dta
* Writes: $figures/maup_sensitivity_<run_tag>.png
* Notes:  Units are ordered by their count. The tract slope is drawn as a dashed reference line because it is the finest unit and the target the aggregations are judged against.

capture confirm file "$tables/maup_${run_tag}.dta"
if _rc {
    display as text "fig_maup_sensitivity: missing maup_${run_tag}.dta; skipping."
    exit 0
}

use "$tables/maup_${run_tag}.dta", clear

gsort -n_units
gen _pos = _n
quietly summarize _pos
local nlev = r(max)

quietly summarize b if key == "tract", meanonly
local b_base = r(mean)

local xlab ""
forvalues i = 1/`nlev' {
    quietly levelsof levlabel if _pos == `i', local(nm) clean
    local xlab `xlab' `i' `"`nm'"'
}

quietly summarize lo
local ymin = r(min)
quietly summarize hi
local ymax = r(max)
local pad = (`ymax' - `ymin') * 0.12
local ylo = `ymin' - `pad'
local yhi = `ymax' + `pad'

twoway (rcap hi lo _pos, lcolor(gs8) lwidth(medthin)) ///
    (scatter b _pos if type == "Baseline", ///
    mcolor(navy)         msymbol(circle)   msize(large)) ///
    (scatter b _pos if type == "Functional", ///
    mcolor(forest_green) msymbol(diamond)  msize(large)) ///
    (scatter b _pos if type == "Administrative", ///
    mcolor(maroon)       msymbol(square)   msize(large)), ///
    yline(`b_base', lpattern(dash) lcolor(navy)) ///
    ytitle("OLS slope: $ mean earnings per +1pp BA share") ///
    xtitle("") ///
    xlabel(`xlab', angle(30) labsize(small) noticks) ///
    xscale(range(0.5 `=`nlev'+0.5')) ///
    yscale(range(`ylo' `yhi')) ///
    ylabel(, format(%9.0gc) angle(0)) ///
    legend(order(2 "Tract (baseline)" 3 "Functional (this study)" ///
    4 "Administrative") position(6) ring(1) cols(3) ///
    region(lstyle(none)) size(small)) ///
    name(_maup, replace)

graph export "$figures/maup_sensitivity_${run_tag}.png", replace width(1600)
display as result "  -> maup_sensitivity_${run_tag}.png"
