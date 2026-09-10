clear all
cap log close

* fig_leiden_sweep_diagnostic.do -- plot the estimated national cluster count and mean modularity against the Leiden CPM resolution, with the chosen gamma and the HC target k marked, one figure per geography.
* Called by: _master.do (generate phase, do_generate switch).
* Reads:  includes _inc/calib_params.do; $tables/leiden_resolution_sweep_<lvl>.csv
* Writes: $figures/leiden_sweep_diagnostic_<lvl>.png
* Notes:  The reference lines come from the constants in _inc/calib_params.do, not from chosen_params_<lvl>.dta, so they can differ from the values calib_04_pick_params.do selected.

include "$program/_inc/calib_params.do"

foreach lvl of global levels {
    capture confirm file "$tables/leiden_resolution_sweep_`lvl'.csv"
    if _rc {
        display as text "  no sweep CSV for `lvl'; skipping."
        continue
    }

    local picked_g = ``lvl'_leiden_res'
    local hc_tgt   = ``lvl'_hc_k_national'

    import delimited using "$tables/leiden_resolution_sweep_`lvl'.csv", varnames(1) clear
    sort resolution

    drop if resolution <= 0

    local title "Leiden CPM resolution diagnostic ({it:`lvl'})"
    local sub ""
    if `picked_g' < . {
        local picked_g_str : display %5.3f `picked_g'
        local hc_tgt_str   : display %12.0fc `hc_tgt'
        local sub "chosen {&gamma} = `=trim("`picked_g_str'")'   (dotted line = HC target of `=trim("`hc_tgt_str'")' clusters)"
    }

    local xl ""
    if `picked_g' < . local xl `"xline(`picked_g', lpattern(dot) lcolor(green))"'

    twoway (line n_clusters_est_national resolution, ///
        sort lcolor(navy) lwidth(medthick) yaxis(1)) ///
        (scatter n_clusters_est_national resolution, ///
        mcolor(navy) msymbol(circle) msize(small) yaxis(1)) ///
        (line mean_modularity resolution, ///
        sort lcolor(maroon) lwidth(medthick) lpattern(dash) yaxis(2)) ///
        (scatter mean_modularity resolution, ///
        mcolor(maroon) msymbol(triangle) msize(small) yaxis(2)), ///
        xtitle("Leiden CPM resolution {&gamma} (log scale)") ///
        ytitle("Estimated national n_clusters", axis(1)) ///
        ytitle("Mean modularity", axis(2)) ///
        title("`title'") subtitle("`sub'") ///
        xscale(log) yscale(log axis(1)) ///
        yline(`hc_tgt', axis(1) lpattern(dot) lcolor(gs8)) ///
        `xl' ///
        legend(order(1 "Estimated n_clusters (left axis)" ///
        3 "Mean modularity (right axis)") ///
        position(6) ring(1) cols(2) region(lstyle(none)))

    graph export "$figures/leiden_sweep_diagnostic_`lvl'.png", replace width(1200)
    display as result "  -> leiden_sweep_diagnostic_`lvl'.png"
}
