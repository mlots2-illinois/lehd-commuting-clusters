clear all
cap log close

* fig_elbow_overlay.do -- overlay the HC and Leiden elbow curves from assess_01_elbow.do on one set of axes, one figure per geography.
* Called by: _master.do (generate phase, do_generate switch).
* Reads:  $temp/elbow_<run_tag>_<lvl>.dta
* Writes: $figures/elbow_overlay_<lvl>.png

foreach lvl of global levels {
    capture confirm file "$temp/elbow_${run_tag}_`lvl'.dta"
    if _rc {
        display as text "  no elbow data for `lvl'; skipping."
        continue
    }
    use "$temp/elbow_${run_tag}_`lvl'.dta", clear
    drop if mi(wcsd_mean)
    sort method n_clusters

    quietly count
    if r(N) < 2 {
        display as text "  elbow_overlay (`lvl'): <2 plottable points; skipping."
        continue
    }

    twoway (line wcsd_mean n_clusters if method == "HC", ///
        sort lcolor(navy) lwidth(medthick)) ///
        (scatter wcsd_mean n_clusters if method == "HC", ///
        mcolor(navy) msymbol(circle) msize(small)) ///
        (line wcsd_mean n_clusters if method == "Leiden", ///
        sort lcolor(maroon) lwidth(medthick)) ///
        (scatter wcsd_mean n_clusters if method == "Leiden", ///
        mcolor(maroon) msymbol(triangle) msize(small)), ///
        xtitle("Number of clusters") ///
        ytitle("Mean within-cluster dissimilarity") ///
        title("Elbow comparison: HC vs Leiden ({it:`lvl'})") ///
        xlabel(, format(%9.0gc)) ///
        legend(order(1 "HC" 3 "Leiden") position(2) ring(0) ///
        cols(1) region(lstyle(none)))
    graph export "$figures/elbow_overlay_`lvl'.png", replace width(1200)
    display as result "  -> elbow_overlay_`lvl'.png"
}
