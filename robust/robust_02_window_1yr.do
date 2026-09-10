* robust_02_window_1yr.do -- rerun the full cluster pipeline on each single year 2019-2023 for the year-to-year stability checks in assess_03_temporal.do.
* Called by: _master.do (robust phase, do_robust_window switch).
* Reads:  includes _inc/robust_prologue.do, _inc/run_cluster.do, _inc/robust_epilogue.do; via run_cluster, $temp/pflow_<lvl>_1yr_<year>.dta and $temp/chosen_params_<lvl>.dta
* Writes: via run_cluster, $temp/final_assignments_temporal_1yr_<year>_<lvl>.dta and the other cluster_* intermediates under the same run_tag
* Notes:  We reuse the parameters calibrated on the 5-year pool for every single year rather than recalibrating, so that the comparison isolates the data window from the parameter choice.

include "$program/_inc/robust_prologue.do"

capture noisily {

    foreach y of numlist 2019/2023 {

        global run_tag             "temporal_1yr_`y'"
        global period_pflow_suffix "1yr_`y'"

        display as result _n "{hline 70}"
        display as result "robust_02: single-year rerun `y' (run_tag=$run_tag, pflow=$period_pflow_suffix)"
        display as result "{hline 70}"

        include "$program/_inc/run_cluster.do"
    }
}
local _body_rc = _rc

include "$program/_inc/robust_epilogue.do"
if `_body_rc' {
    display as error "robust_02: body failed (rc=`_body_rc'); globals restored."
    exit `_body_rc'
}
