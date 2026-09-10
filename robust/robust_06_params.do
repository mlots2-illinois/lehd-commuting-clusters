* robust_06_params.do -- rerun the cluster pipeline with each calibrated parameter (alpha, cap_km, hc_k_national, leiden_resolution) scaled by 0.5, 0.75, 1.25 and 1.5 in turn.
* Called by: _master.do (robust phase, do_robust_params switch).
* Reads:  includes _inc/robust_prologue.do, _inc/sens_grid.do, _inc/run_cluster.do, _inc/robust_epilogue.do; $temp/chosen_params_<lvl>.dta
* Writes: $temp/chosen_params_<lvl>.dta (perturbed, then restored); $temp/_chosen_params_baseline_<lvl>.dta (backup, erased at the end); via run_cluster, $temp/final_assignments_sens_<p>_m<mult>_<lvl>.dta
* Notes:  The cluster scripts read parameters only from chosen_params_<lvl>.dta, so we overwrite that file for each specification and restore it from a backup afterwards. A stale backup at start means the previous run crashed; we restore from it before doing anything else.

include "$program/_inc/robust_prologue.do"

foreach lvl of global levels {
    capture confirm file "$temp/_chosen_params_baseline_`lvl'.dta"
    if !_rc {
        display as text "robust_06: WARNING: stale backup _chosen_params_baseline_`lvl'.dta found" ///
            " (previous run crashed?); restoring chosen_params_`lvl'.dta from it."
        copy "$temp/_chosen_params_baseline_`lvl'.dta" "$temp/chosen_params_`lvl'.dta", replace
    }
}

foreach lvl of global levels {
    copy "$temp/chosen_params_`lvl'.dta" "$temp/_chosen_params_baseline_`lvl'.dta", replace
}

include "$program/_inc/sens_grid.do"

local n_specs : word count `sens_tags'

capture noisily {

    forvalues i = 1/`n_specs' {
        local rt    : word `i' of `sens_tags'
        local pname : word `i' of `sens_params'
        local v_t   : word `i' of `sens_v_t'
        local v_c   : word `i' of `sens_v_c'

        display as result _n "{hline 70}"
        display as result "robust_06: `rt'  (`pname': tract=`v_t' county=`v_c')"
        display as result "{hline 70}"

        foreach lvl of global levels {
            local vv = cond("`lvl'" == "tract", "`v_t'", "`v_c'")
            use "$temp/_chosen_params_baseline_`lvl'.dta", clear
            replace `pname' = `vv'
            save "$temp/chosen_params_`lvl'.dta", replace
        }

        global run_tag "`rt'"
        include "$program/_inc/run_cluster.do"
    }
}
local _body_rc = _rc

foreach lvl of global levels {
    capture confirm file "$temp/_chosen_params_baseline_`lvl'.dta"
    if !_rc {
        copy  "$temp/_chosen_params_baseline_`lvl'.dta" "$temp/chosen_params_`lvl'.dta", replace
        erase "$temp/_chosen_params_baseline_`lvl'.dta"
    }
}

include "$program/_inc/robust_epilogue.do"
if `_body_rc' {
    display as error "robust_06: body failed (rc=`_body_rc'); chosen_params and globals restored."
    exit `_body_rc'
}
