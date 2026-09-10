clear all
cap log close

* fig_maup_paramsweep.do -- plot the earnings-on-BA slope across each parameter perturbation from assess_10_maup_paramsweep.do, four panels with the tract and county slopes as references.
* Called by: _master.do (generate phase, do_generate switch).
* Reads:  $tables/maup_paramsweep_<run_tag>.dta
* Writes: $figures/maup_paramsweep_<run_tag>.png
* Notes:  The panel titles hard-code the calibrated tract values (0.80, 100 km, 9,000, 0.171). grc1leg is installed from SSC if absent so the four panels share one legend; graph combine is the fallback.

capture confirm file "$tables/maup_paramsweep_${run_tag}.dta"
if _rc {
    display as text "fig_maup_paramsweep: missing maup_paramsweep_${run_tag}.dta; skipping."
    exit 0
}

use "$tables/maup_paramsweep_${run_tag}.dta", clear

quietly summarize lo
local ymin = r(min)
quietly summarize hi
local ymax = r(max)
quietly summarize ref_county, meanonly
local ymin = min(`ymin', r(mean))
quietly summarize ref_tract, meanonly
local ymax = max(`ymax', r(mean))
local pad = (`ymax' - `ymin') * 0.08
local ylo = `ymin' - `pad'
local yhi = `ymax' + `pad'

quietly summarize ref_tract,  meanonly
local rt = r(mean)
quietly summarize ref_county, meanonly
local rc = r(mean)

capture which grc1leg
if _rc capture ssc install grc1leg, replace
capture which grc1leg
local have_grc = (_rc == 0)

local panels ""
foreach param in alpha cap k gamma {
    if "`param'" == "alpha" local ptitle "{&alpha} (chosen 0.80)"
    if "`param'" == "cap"   local ptitle "cap (chosen 100 km)"
    if "`param'" == "k"     local ptitle "{it:k}{sub:HC} (chosen 9,000)"
    if "`param'" == "gamma" local ptitle "{&gamma} (chosen 0.171)"

    local lyr ""
    local lyr `lyr' (rcap hi lo mult if param=="`param'", lcolor(gs10) lwidth(thin))
    quietly count if param=="`param'" & method=="hc"
    if r(N) > 0 ///
        local lyr `lyr' (connected b mult if param=="`param'" & method=="hc", ///
        sort lcolor(navy) mcolor(navy) msymbol(circle) lwidth(medthin))
    quietly count if param=="`param'" & method=="leiden"
    if r(N) > 0 ///
        local lyr `lyr' (connected b mult if param=="`param'" & method=="leiden", ///
        sort lcolor(forest_green) mcolor(forest_green) msymbol(diamond) lwidth(medthin))

    local common ///
        yline(`rt', lpattern(dash) lcolor(gs6)) ///
        yline(`rc', lpattern(shortdash) lcolor(cranberry)) ///
        xline(1, lpattern(dot) lcolor(gs12)) ///
        title("`ptitle'", size(medsmall)) ///
        ytitle("") xtitle("") ///
        xlabel(0.5 0.75 1 1.25 1.5, labsize(small)) ///
        ylabel(800(200)1400, format(%9.0gc) labsize(small) angle(0)) ///
        yscale(range(`ylo' `yhi'))

    if "`param'" == "alpha" {
        twoway `lyr', `common' name(g_`param', replace) nodraw ///
            legend(order(2 "HC" 3 "Leiden") rows(1) region(lstyle(none)) size(small))
    }
    else {
        twoway `lyr', `common' legend(off) name(g_`param', replace) nodraw
    }
    local panels `panels' g_`param'
}

local ttls l1title("OLS slope: $ mean earnings per +1pp BA share", size(small)) ///
    b1title("Multiple of chosen parameter value", size(small))
if `have_grc' {
    grc1leg `panels', cols(2) `ttls' name(g_combined, replace)
}
else {
    graph combine `panels', cols(2) `ttls' name(g_combined, replace)
}

graph export "$figures/maup_paramsweep_${run_tag}.png", replace width(1800)
display as result "  -> maup_paramsweep_${run_tag}.png"
display as text "Note: dashed = tract (finest) slope; short-dash = county (administrative) slope."
