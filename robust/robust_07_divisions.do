* robust_07_divisions.do -- rerun the tract cluster pipeline under six alternative division geometries (Wwide, Wnarrow, ovlarge, ovsmall, bands20, bands45) built in build_03_divisions.do.
* Called by: _master.do (robust phase, do_robust_divisions switch).
* Reads:  includes _inc/robust_prologue.do, _inc/run_cluster.do, _inc/robust_epilogue.do; via run_cluster, $temp/division_lookup_<div_tag>.dta, $temp/pflow_tract_5yr_2019_2023.dta and $temp/chosen_params_tract.dta
* Writes: via run_cluster, $temp/final_assignments_div_<div_tag>_tract.dta and the other cluster_* intermediates under the same run_tag
* Notes:  Only the tract geography is divided into overlapping bands, so we set $levels to "tract" for the body; the epilogue restores it. Six of the fifteen geometries in _inc/div_specs.do are rerun here.

if strpos(" $levels ", " tract ") == 0 {
    display as text "robust_07: 'tract' not in \$levels; the division check is tract-only, skipping."
    exit 0
}

include "$program/_inc/robust_prologue.do"

global levels "tract"

capture noisily {

    foreach g in Wwide Wnarrow ovlarge ovsmall bands20 bands45 {

        global div_tag "`g'"
        global run_tag "div_`g'"

        display as result _n "{hline 70}"
        display as result "robust_07: division geometry `g' (run_tag=$run_tag)"
        display as result "{hline 70}"

        include "$program/_inc/run_cluster.do"
    }
}
local _body_rc = _rc

include "$program/_inc/robust_epilogue.do"
if `_body_rc' {
    display as error "robust_07: body failed (rc=`_body_rc'); globals restored."
    exit `_body_rc'
}
