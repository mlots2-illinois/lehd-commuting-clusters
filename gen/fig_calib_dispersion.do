clear all

* fig_calib_dispersion.do -- draw the dispersion-of-D lines from the alpha and cap sweep of calib_01_alpha_cap.do, one panel per geography.
* Called by: _master.do (generate phase, do_generate switch).
* Reads:  $tables/calib_dispersion_<lvl>.csv; shell "$py" _py/calib_dispersion_lines.py
* Writes: $figures/calib_dispersion_lines_<lvl>.png
* Notes:  Both geographies must be present; the script skips rather than draw one panel, because the paper shows the two side by side.

local pyscript "$program/_py/calib_dispersion_lines.py"

local _missing 0
foreach lvl in tract county {
    capture confirm file "$tables/calib_dispersion_`lvl'.csv"
    if _rc local _missing 1
}
if `_missing' {
    display as text "fig_calib_dispersion: calib dispersion CSV(s) absent (run calib_01 for both levels); skipping."
    exit 0
}

capture erase "$figures/calib_dispersion_lines_tract.png"
capture erase "$figures/calib_dispersion_lines_county.png"
shell "$py" "`pyscript'" "$tables" "$figures"
local _failed 0
foreach lvl in tract county {
    capture confirm file "$figures/calib_dispersion_lines_`lvl'.png"
    if _rc local _failed 1
}
if `_failed' {
    display as error "fig_calib_dispersion: python plotter failed; PNG(s) not produced."
    exit 0
}
display as result "  -> calib_dispersion_lines_{tract,county}.png"
