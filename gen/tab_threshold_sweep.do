clear all

* tab_threshold_sweep.do -- write the minimum-flow threshold table comparing the robust_05_threshold.do reruns (min_flow = 2, 3, 5, 10, 20) with the baseline.
* Called by: _master.do (generate phase, do_generate switch).
* Reads:  includes _inc/booktabs.do and _inc/compare_to_baseline.do; $temp/final_assignments_baseline_<lvl>.dta; $temp/final_assignments_thresh_<m>_<lvl>.dta
* Writes: $tables/threshold_sweep.tex; $tables/threshold_sweep.csv

include "$program/_inc/booktabs.do"
include "$program/_inc/compare_to_baseline.do"

tempname pf
tempfile res
postfile `pf' str6 geo int minflow str8 method double NMI double ARI using "`res'", replace

foreach lvl of global levels {

    capture confirm file "$temp/final_assignments_baseline_`lvl'.dta"
    if _rc continue

    post `pf' ("`lvl'") (0) ("hc")     (1) (1)
    post `pf' ("`lvl'") (0) ("leiden") (1) (1)

    foreach m in 2 3 5 10 20 {
        capture confirm file "$temp/final_assignments_thresh_`m'_`lvl'.dta"
        if _rc {
            display as text "  no thresh_`m' rerun at `lvl'; skipping that row."
            continue
        }
        compare_to_baseline, lvl(`lvl') alt("$temp/final_assignments_thresh_`m'_`lvl'.dta")
        local nh = r(NMI_hc)
        local ah = r(ARI_hc)
        local nl = r(NMI_le)
        local al = r(ARI_le)
        post `pf' ("`lvl'") (`m') ("hc")     (`nh') (`ah')
        post `pf' ("`lvl'") (`m') ("leiden") (`nl') (`al')
    }
}
postclose `pf'

use "`res'", clear
quietly count
if r(N) == 0 {
    display as text "tab_threshold_sweep: nothing to write; skipping."
    exit 0
}
reshape wide NMI ARI, i(geo minflow) j(method) string
gen byte _go = cond(geo == "tract", 1, 2)
gsort _go minflow

export delimited geo minflow NMIhc ARIhc NMIleiden ARIleiden ///
    using "$tables/threshold_sweep.csv", replace

bt_open, handle(_t) path("$tables/threshold_sweep.tex") colspec(ll cc cc) ///
    script(gen/tab_threshold_sweep.do)
file write _t "& & \multicolumn{2}{c}{\textbf{Hierarchical}} & \multicolumn{2}{c}{\textbf{Leiden}} \\" _n
file write _t "\cmidrule(lr){3-4} \cmidrule(lr){5-6}" _n
file write _t "\textbf{Geography} & \textbf{Min.\ flow} & NMI & ARI & NMI & ARI \\" _n
file write _t "\midrule" _n
local prev ""
forvalues i = 1/`=_N' {
    local g    = geo[`i']
    local glab = cond("`g'" != "`prev'", strproper("`g'"), "")
    if "`g'" != "`prev'" & "`prev'" != "" file write _t "\addlinespace" _n
    local prev "`g'"
    local mf = minflow[`i']
    if `mf' == 0 local mfc "\textbf{0}\,$^\dagger$"
    else         local mfc "`mf'"
    bt_num, value(`=NMIhc[`i']') fmt(%5.3f)
    local nh `"`r(s)'"'
    bt_num, value(`=ARIhc[`i']') fmt(%5.3f)
    local ah `"`r(s)'"'
    bt_num, value(`=NMIleiden[`i']') fmt(%5.3f)
    local nl `"`r(s)'"'
    bt_num, value(`=ARIleiden[`i']') fmt(%5.3f)
    local al `"`r(s)'"'
    bt_row, handle(_t) cells(`" "`glab'" "`mfc'" "`nh'" "`ah'" "`nl'" "`al'" "')
}
bt_close, handle(_t)

display as result "  -> threshold_sweep.tex + threshold_sweep.csv"
