clear all
cap log close

* fig_elbow_hc_clusters.do -- plot mean within-cluster dissimilarity against the number of HC clusters from the assess_01_elbow.do grid, one figure per geography.
* Called by: _master.do (generate phase, do_generate switch).
* Reads:  $temp/elbow_<run_tag>_<lvl>.dta
* Writes: $figures/elbow_hc_clusters_vs_wcsd_<lvl>.png

foreach lvl of global levels {
    capture confirm file "$temp/elbow_${run_tag}_`lvl'.dta"
    if _rc {
        display as text "  no elbow data for `lvl'; skipping."
        continue
    }
    use "$temp/elbow_${run_tag}_`lvl'.dta", clear
    keep if method == "HC"
    drop if mi(wcsd_mean)
    sort n_clusters

    quietly count
    if r(N) < 2 {
        display as text "  elbow_hc (`lvl'): <2 plottable points; skipping."
        continue
    }

    twoway (line wcsd_mean n_clusters, sort lcolor(black) mcolor(black)) ///
        (scatter wcsd_mean n_clusters, mcolor(black) msymbol(circle)), ///
        xtitle("Number of clusters") ///
        ytitle("Mean within-cluster dissimilarity") ///
        title("HC elbow: n_clusters vs mean within-cluster D (`lvl')") ///
        xlabel(, format(%9.0gc)) ///
        legend(off)
    graph export "$figures/elbow_hc_clusters_vs_wcsd_`lvl'.png", replace width(1200)
    display as result "  -> elbow_hc_clusters_vs_wcsd_`lvl'.png"
}
