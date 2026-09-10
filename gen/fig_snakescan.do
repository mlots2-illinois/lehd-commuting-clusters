clear all

* fig_snakescan.do -- draw the appendix figure illustrating the snake-scan division geometry.
* Called by: _master.do (generate phase, do_generate switch).
* Reads:  $tables/leiden_resolution_sweep_tract.csv; shell "$py" _py/appendix_figures.py ... snakescan
* Writes: $figures/appendix_snakescan.png

local pyscript "$program/_py/appendix_figures.py"
local sweepcsv "$tables/leiden_resolution_sweep_tract.csv"

capture erase "$figures/appendix_snakescan.png"
shell "$py" "`pyscript'" "`sweepcsv'" "$figures" snakescan
capture confirm file "$figures/appendix_snakescan.png"
if _rc {
    display as error "fig_snakescan: python plotter failed; no PNG produced."
    exit 0
}
display as result "  -> appendix_snakescan.png"
