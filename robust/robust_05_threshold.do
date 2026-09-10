* robust_05_threshold.do -- rerun the cluster pipeline on the minimum-flow pooled flows (min_flow = 2, 3, 5, 10, 20) built in build_01_pflows.do.
* Called by: _master.do (robust phase, do_robust_threshold switch).
* Reads:  includes _inc/robust_prologue.do, _inc/run_cluster.do, _inc/robust_epilogue.do; via run_cluster, $temp/pflow_<lvl>_5yr_2019_2023_mf<m>.dta and $temp/chosen_params_<lvl>.dta
* Writes: via run_cluster, $temp/final_assignments_thresh_<m>_<lvl>.dta and the other cluster_* intermediates under the same run_tag
* Notes:  The threshold changes only the edge list; the calibrated parameters stay fixed so the comparison in tab_threshold_sweep.do reflects the dropped small flows alone.

include "$program/_inc/robust_prologue.do"

capture noisily {

    foreach m in 2 3 5 10 20 {

        global min_flow            `m'
        global flow_tag            "_mf`m'"
        global run_tag             "thresh_`m'"
        global period_pflow_suffix "5yr_${yr_start}_${yr_end}_mf`m'"

        display as result _n "{hline 70}"
        display as result "robust_05: minimum flow >= `m' (run_tag=$run_tag)"
        display as result "{hline 70}"

        include "$program/_inc/run_cluster.do"
    }
}
local _body_rc = _rc

include "$program/_inc/robust_epilogue.do"
if `_body_rc' {
    display as error "robust_05: body failed (rc=`_body_rc'); globals restored."
    exit `_body_rc'
}
