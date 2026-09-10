clear all
cap log close

* fig_sensitivity_dotplot.do -- dot plot of ARI against the baseline for every parameter perturbation in assess_04_sensitivity.do, one figure per geography.
* Called by: _master.do (generate phase, do_generate switch).
* Reads:  $tables/sensitivity_summary.dta
* Writes: $figures/sensitivity_dotplot_<lvl>.png
* Notes:  Rows are grouped by parameter (alpha, cap, k, gamma) with a half-row gap between groups so the four sweeps read as blocks.

capture confirm file "$tables/sensitivity_summary.dta"
if _rc {
    display as text "fig_sensitivity_dotplot: missing sensitivity_summary.dta; skipping."
    exit 0
}

foreach lvl of global levels {
    use "$tables/sensitivity_summary.dta", clear
    keep if geo_level == "`lvl'"
    quietly count
    if r(N) == 0 continue

    gen str30 perturb = ""
    replace perturb = "{&alpha} = " + string(value, "%4.2f")         if param == "alpha"
    replace perturb = "cap = "      + string(value, "%5.0f") + " km"  if param == "cap_km"
    replace perturb = "{it:k} = "   + string(value, "%6.0fc")         if param == "hc_k_national"
    replace perturb = "{&gamma} = " + string(value, "%5.3f")         if param == "leiden_resolution"

    gen byte _prank = .
    replace _prank = 1 if param == "alpha"
    replace _prank = 2 if param == "cap_km"
    replace _prank = 3 if param == "hc_k_national"
    replace _prank = 4 if param == "leiden_resolution"

    egen _order  = group(_prank value)
    egen _grpidx = group(_prank)
    quietly summarize _order
    local nrows = r(max)
    gen double _ypos = _order + 0.5 * (_grpidx - 1)

    local ylab ""
    forvalues i = 1/`nrows' {
        quietly summarize _ypos if _order == `i', meanonly
        local p = r(mean)
        quietly levelsof perturb if _order == `i', local(plab) clean
        if `"`plab'"' != "" local ylab `ylab' `p' `"`plab'"'
    }
    quietly summarize _ypos
    local ybot = r(min) - 1
    local ytop = r(max) + 1

    twoway (scatter _ypos ARI if method == "hc", ///
        mcolor(navy) msymbol(circle) msize(medlarge)) ///
        (scatter _ypos ARI if method == "leiden", ///
        mcolor(maroon) msymbol(triangle) msize(medlarge)), ///
        xtitle("ARI vs. baseline partition") ytitle("") ///
        title("Sensitivity to calibration perturbations ({it:`lvl'})") ///
        xline(1, lpattern(dot) lcolor(gs10)) xlabel(0(0.2)1.0) ///
        ylabel(`ylab', angle(0) labsize(small) noticks) ///
        yscale(range(`ybot' `ytop')) ///
        legend(order(1 "HC" 2 "Leiden") position(6) ring(1) ///
        cols(2) region(lstyle(none)))

    graph export "$figures/sensitivity_dotplot_`lvl'.png", replace width(1200)
    display as result "  -> sensitivity_dotplot_`lvl'.png"
}
