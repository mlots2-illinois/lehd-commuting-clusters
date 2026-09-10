* robust_03_window_baseline_ref.do -- placeholder for the 5-year window, which is the baseline run itself; it prints a notice and reruns nothing.
* Called by: _master.do (robust phase, do_robust_window switch).
* Reads:  includes _inc/robust_prologue.do and _inc/robust_epilogue.do
* Notes:  The file exists so the window series 01-03 is complete in the timing log; final_assignments_baseline_<lvl>.dta already serves as the 5-year result.

include "$program/_inc/robust_prologue.do"

capture noisily {

    display as result _n "{hline 70}"
    display as result "robust_03: 5-year window = baseline (final_assignments_baseline_*); nothing to re-run."
    display as result "{hline 70}"
}
local _body_rc = _rc

include "$program/_inc/robust_epilogue.do"
if `_body_rc' {
    display as error "robust_03: body failed (rc=`_body_rc'); globals restored."
    exit `_body_rc'
}
