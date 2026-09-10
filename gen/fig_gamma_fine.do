clear all

* fig_gamma_fine.do -- draw the fine Leiden resolution sweep from calib_02_resolution.do with the chosen gamma marked, one figure per geography.
* Called by: _master.do (generate phase, do_generate switch).
* Reads:  includes _inc/get_chosen.do; $tables/leiden_resolution_sweep_<lvl>.csv; $temp/chosen_params_<lvl>.dta; shell "$py" _py/leiden_resolution_lines.py
* Writes: $figures/leiden_resolution_fine_<lvl>.png
* Notes:  The plotter may still write leiden_resolution_lines_<lvl>.png under its older name; in that case we copy it to the manuscript name and erase the original.

include "$program/_inc/get_chosen.do"

local pyscript "$program/_py/leiden_resolution_lines.py"
foreach lvl of global levels {
    capture confirm file "$tables/leiden_resolution_sweep_`lvl'.csv"
    if _rc {
        display as text "fig_gamma_fine: no sweep CSV for `lvl'; skipping."
        continue
    }

    local pickopt ""
    get_chosen, lvl(`lvl') param(leiden_resolution)
    if r(value) < . local pickopt `"--picked-gamma `=r(value)'"'

    capture erase "$figures/leiden_resolution_fine_`lvl'.png"
    capture erase "$figures/leiden_resolution_lines_`lvl'.png"
    shell "$py" "`pyscript'" "$tables" "$figures" "`lvl'" `pickopt'

    capture confirm file "$figures/leiden_resolution_fine_`lvl'.png"
    if _rc {
        capture confirm file "$figures/leiden_resolution_lines_`lvl'.png"
        if !_rc {
            display as text "  note: plotter wrote leiden_resolution_lines_`lvl'.png; copying to the manuscript name."
            copy "$figures/leiden_resolution_lines_`lvl'.png" ///
                "$figures/leiden_resolution_fine_`lvl'.png", replace
        }
    }

    capture confirm file "$figures/leiden_resolution_fine_`lvl'.png"
    if _rc {
        display as error "fig_gamma_fine: python plotter failed for `lvl'; no PNG produced."
        continue
    }
    capture erase "$figures/leiden_resolution_lines_`lvl'.png"
    display as result "  -> leiden_resolution_fine_`lvl'.png"
}
