* robust_08_alpha_full.do -- rerun the cluster pipeline at alpha = 0 and alpha = 0.2, the endpoints the multiplier grid in robust_06_params.do cannot reach.
* Called by: _master.do (robust phase, do_robust_alpha_full switch).
* Reads:  includes _inc/robust_prologue.do, _inc/run_cluster.do, _inc/robust_epilogue.do; $temp/chosen_params_<lvl>.dta
* Writes: $temp/chosen_params_<lvl>.dta (perturbed, then restored); $temp/_chosen_params_baseline_<lvl>.dta (backup, erased at the end); via run_cluster, $temp/final_assignments_sens_a_full<aaa>_<lvl>.dta
* Notes:  The grid is $alpha_full_grid (default "0.0 0.2") and the tag is alpha x 100, zero-padded to three digits. Together with the robust_06 runs and the baseline this gives alpha on a 0 to 1 grid in steps of 0.2.
*
* robust_06_params.do perturbs each calibrated parameter by a multiplier, half to
* one-and-a-half times its value. That is the natural perturbation for the cap,
* k and gamma, which have no intrinsic upper bound. Alpha is a weight on [0, 1],
*     D_ij = alpha (1 - P_ij) + (1 - alpha) min(d_ij / dbar, 1),
* so a multiplier misrepresents it and cannot reach its endpoints. At alpha = 0
* the partition is built from distance alone; at alpha = 1 distance drops out.
* This script supplies the two settings the multiplier grid cannot, so that
* alpha can be reported on a 0 to 1 grid in steps of 0.2:
*     alpha 0.0  -> sens_a_full000   (this script)
*     alpha 0.2  -> sens_a_full020   (this script)
*     alpha 0.4  -> sens_a_m050      (robust_06, 0.5 x 0.8)
*     alpha 0.6  -> sens_a_m075      (robust_06, 0.75 x 0.8)
*     alpha 0.8  -> baseline         (the calibrated value)
*     alpha 1.0  -> sens_a_m125      (robust_06, 1.25 x 0.8)
* The script loops over $levels, so both geographies run unless $levels is
* narrowed before it is called.

include "$program/_inc/robust_prologue.do"

if "$alpha_full_grid" == "" global alpha_full_grid "0.0 0.2"

* Restore any perturbed parameter file left behind by a crashed run before
* taking a fresh backup, exactly as robust_06_params.do does.
foreach lvl of global levels {
    capture confirm file "$temp/_chosen_params_baseline_`lvl'.dta"
    if !_rc {
        display as text "robust_08: WARNING: stale backup for `lvl' found" ///
            " (previous run crashed?); restoring chosen_params_`lvl'.dta."
        copy "$temp/_chosen_params_baseline_`lvl'.dta" "$temp/chosen_params_`lvl'.dta", replace
    }
}
foreach lvl of global levels {
    copy "$temp/chosen_params_`lvl'.dta" "$temp/_chosen_params_baseline_`lvl'.dta", replace
}

capture noisily {
    foreach a of global alpha_full_grid {
        local atag = string(`a' * 100, "%03.0f")
        local rt "sens_a_full`atag'"

        display as result _n "{hline 70}"
        display as result "robust_08: `rt'  (alpha = `a')"
        display as result "{hline 70}"

        foreach lvl of global levels {
            use "$temp/_chosen_params_baseline_`lvl'.dta", clear
            replace alpha = `a'
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
    display as error "robust_08: run failed with rc=`_body_rc'; chosen_params and globals restored."
    exit `_body_rc'
}
