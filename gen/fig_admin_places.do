clear all

* fig_admin_places.do -- draw the appendix figure comparing the Leiden tract sweep against administrative place counts.
* Called by: _master.do (generate phase, do_generate switch).
* Reads:  $tables/leiden_resolution_sweep_tract.csv; shell "$py" _py/appendix_figures.py ... admin
* Writes: $figures/appendix_admin_vs_places.png

local pyscript "$program/_py/appendix_figures.py"
local sweepcsv "$tables/leiden_resolution_sweep_tract.csv"

capture erase "$figures/appendix_admin_vs_places.png"
shell "$py" "`pyscript'" "`sweepcsv'" "$figures" admin
capture confirm file "$figures/appendix_admin_vs_places.png"
if _rc {
    display as error "fig_admin_places: python plotter failed; no PNG produced."
    exit 0
}
display as result "  -> appendix_admin_vs_places.png"
