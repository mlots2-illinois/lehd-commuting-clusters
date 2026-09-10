* robust_01_window_3yr.do -- rerun the full cluster pipeline on the 2021-2023 pooled flows so the assess phase can compare a 3-year window with the 5-year baseline.
* Called by: _master.do (robust phase, do_robust_window switch).
* Reads:  includes _inc/robust_prologue.do, _inc/run_cluster.do, _inc/robust_epilogue.do; via run_cluster, $temp/pflow_<lvl>_3yr_2021_2023.dta and $temp/chosen_params_<lvl>.dta
* Writes: via run_cluster, $temp/final_assignments_temporal_3yr_2021_2023_<lvl>.dta and the other cluster_* intermediates under the same run_tag
* Notes:  The prologue stashes the baseline globals and the epilogue restores them, so a failed body still leaves _master.do in the baseline namespace; the script then exits with the body's return code.

include "$program/_inc/robust_prologue.do"

global run_tag             "temporal_3yr_2021_2023"
global period_pflow_suffix "3yr_2021_2023"

capture noisily {

    display as result _n "{hline 70}"
    display as result "robust_01: 3-year window rerun (run_tag=$run_tag, pflow=$period_pflow_suffix)"
    display as result "{hline 70}"

    include "$program/_inc/run_cluster.do"
}
local _body_rc = _rc

include "$program/_inc/robust_epilogue.do"
if `_body_rc' {
    display as error "robust_01: body failed (rc=`_body_rc'); globals restored."
    exit `_body_rc'
}
