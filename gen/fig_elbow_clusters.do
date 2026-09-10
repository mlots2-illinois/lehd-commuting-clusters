clear all
cap log close

* fig_elbow_clusters.do -- draw the HC and Leiden elbow curves from assess_01_elbow.do as two side-by-side panels, one figure per geography.
* Called by: _master.do (generate phase, do_generate switch).
* Reads:  $temp/elbow_<run_tag>_<lvl>.dta
* Writes: $figures/elbow_clusters_vs_wcsd_<lvl>.png
* Notes:  The combined figure is written only when both methods have at least two plottable points; a one-sided figure would invite a comparison the data cannot support.

foreach lvl of global levels {
    capture confirm file "$temp/elbow_${run_tag}_`lvl'.dta"
    if _rc {
        display as text "  no elbow data for `lvl'; skipping."
        continue
    }
    use "$temp/elbow_${run_tag}_`lvl'.dta", clear
    drop if mi(wcsd_mean)

    local _haveHC 0
    preserve
    keep if method == "HC"
    sort n_clusters
    quietly count
    if r(N) >= 2 {
        local _haveHC 1
        twoway (line wcsd_mean n_clusters, sort lcolor(navy) lwidth(medthick)) ///
            (scatter wcsd_mean n_clusters, mcolor(navy) msymbol(circle) msize(small)), ///
            xtitle("Number of clusters") ///
            ytitle("Mean within-cluster dissimilarity") ///
            title("Hierarchical clustering") ///
            xlabel(, format(%9.0gc)) legend(off) name(_elbHC, replace)
    }
    restore

    local _haveLE 0
    preserve
    keep if method == "Leiden"
    sort n_clusters
    quietly count
    if r(N) >= 2 {
        local _haveLE 1
        twoway (line wcsd_mean n_clusters, sort lcolor(maroon) lwidth(medthick)) ///
            (scatter wcsd_mean n_clusters, mcolor(maroon) msymbol(triangle) msize(small)), ///
            xtitle("Number of clusters") ///
            ytitle("Mean within-cluster dissimilarity") ///
            title("Leiden community detection") ///
            xlabel(, format(%9.0gc)) legend(off) name(_elbLE, replace)
    }
    restore

    if `_haveHC' & `_haveLE' {
        graph combine _elbHC _elbLE, rows(1) ///
            title("Within-cluster dissimilarity vs number of clusters ({it:`lvl'})") ///
            name(_elbBOTH, replace)
        graph export "$figures/elbow_clusters_vs_wcsd_`lvl'.png", replace width(2000)
        display as result "  -> elbow_clusters_vs_wcsd_`lvl'.png"
    }
    else {
        display as text "  elbow_clusters (`lvl'): a method had <2 plottable points; skipping combine."
    }
    capture graph drop _elbHC _elbLE _elbBOTH
}
