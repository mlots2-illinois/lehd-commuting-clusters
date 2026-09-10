* get_chosen.do -- defines get_chosen, which returns one chosen calibration parameter for a level from the persisted chosen_params file, or missing.
* Included by: gen/fig_gamma_fine.do, tab_chosen_params.do, tab_leiden_sweep.do, tab_dbar_sweep.do, tab_alpha_sweep.do.
* Expects: $temp
* Reads:  $temp/chosen_params_<lvl>.dta (written by calib/calib_04_pick_params.do)
* Notes:  the file is read in a temporary frame so the caller's data are untouched.

capture program drop get_chosen
program define get_chosen, rclass
    syntax , lvl(string) param(name)

    local v = .
    capture confirm file "$temp/chosen_params_`lvl'.dta"
    if !_rc {
        tempname fr
        frame create `fr'
        frame `fr' {
            use "$temp/chosen_params_`lvl'.dta", clear
            capture confirm variable `param'
            if !_rc local v = `param'[1]
        }
        frame drop `fr'
    }
    return scalar value = `v'
end
