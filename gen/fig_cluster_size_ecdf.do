clear all
cap log close

* fig_cluster_size_ecdf.do -- plot the empirical CDF of cluster size (units per cluster) for HC and Leiden on a log axis, one figure per geography.
* Called by: _master.do (generate phase, do_generate switch).
* Reads:  includes _inc/cluster_sizes.do; $temp/final_assignments_<run_tag>_<lvl>.dta
* Writes: $figures/cluster_size_ecdf_<lvl>.png

include "$program/_inc/cluster_sizes.do"

foreach lvl of global levels {
    capture confirm file "$temp/final_assignments_${run_tag}_`lvl'.dta"
    if _rc {
        display as text "  no final_assignments for `lvl'; skipping."
        continue
    }
    use geoid hc_cluster leiden_cluster using ///
        "$temp/final_assignments_${run_tag}_`lvl'.dta", clear

    cluster_sizes

    sort method size
    bysort method: gen long _rank = _n
    bysort method: gen long _n_total = _N
    gen double ecdf = _rank / _n_total

    twoway (line ecdf size if method == "HC", ///
        sort lcolor(navy) lwidth(medthick)) ///
        (line ecdf size if method == "Leiden", ///
        sort lcolor(maroon) lwidth(medthick) lpattern(dash)), ///
        xtitle("Cluster size (units per cluster, log scale)") ///
        ytitle("Empirical CDF") ///
        title("Cluster size distribution ({it:`lvl'})") ///
        xscale(log) ///
        xlabel(1 5 25 100 500 2000, format(%9.0gc)) ///
        legend(order(1 "HC" 2 "Leiden") position(5) ring(0) ///
        cols(1) region(lstyle(none)))

    graph export "$figures/cluster_size_ecdf_`lvl'.png", replace width(1200)
    display as result "  -> cluster_size_ecdf_`lvl'.png"
}
