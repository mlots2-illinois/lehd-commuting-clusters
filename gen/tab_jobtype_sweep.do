clear all

* tab_jobtype_sweep.do -- write the job-type table comparing the primary-jobs (JT01) rerun with the all-jobs (JT00) baseline.
* Called by: _master.do (generate phase, do_generate switch).
* Reads:  includes _inc/booktabs.do and _inc/compare_to_baseline.do; $temp/final_assignments_baseline_<lvl>.dta; $temp/final_assignments_jt01_<lvl>.dta
* Writes: $tables/jobtype_sweep.tex; $tables/jobtype_sweep.csv
* Notes:  The JT01 rerun comes from the do_robust_jobtype block of _master.do, which calls _inc/run_cluster.do directly rather than a robust_*.do script.

include "$program/_inc/booktabs.do"
include "$program/_inc/compare_to_baseline.do"

tempname pf
tempfile res
postfile `pf' str6 geo str8 jobtype str8 method double NMI double ARI using "`res'", replace

foreach lvl of global levels {

    capture confirm file "$temp/final_assignments_baseline_`lvl'.dta"
    if _rc {
        display as text "  no baseline at `lvl'; skipping."
        continue
    }

    post `pf' ("`lvl'") ("JT00") ("hc")     (1) (1)
    post `pf' ("`lvl'") ("JT00") ("leiden") (1) (1)

    capture confirm file "$temp/final_assignments_jt01_`lvl'.dta"
    if _rc {
        display as text "  no JT01 run at `lvl' (build JT01 + run cluster phase); baseline row only."
        continue
    }

    compare_to_baseline, lvl(`lvl') alt("$temp/final_assignments_jt01_`lvl'.dta")
    local nh = r(NMI_hc)
    local ah = r(ARI_hc)
    local nl = r(NMI_le)
    local al = r(ARI_le)
    post `pf' ("`lvl'") ("JT01") ("hc")     (`nh') (`ah')
    post `pf' ("`lvl'") ("JT01") ("leiden") (`nl') (`al')
}
postclose `pf'

use "`res'", clear
quietly count
if r(N) == 0 {
    display as text "tab_jobtype_sweep: nothing to write; skipping."
    exit 0
}
reshape wide NMI ARI, i(geo jobtype) j(method) string
gen byte _go  = cond(geo == "tract", 1, 2)
gen byte _jo  = cond(jobtype == "JT00", 1, 2)
gsort _go _jo

export delimited geo jobtype NMIhc ARIhc NMIleiden ARIleiden ///
    using "$tables/jobtype_sweep.csv", replace

bt_open, handle(_t) path("$tables/jobtype_sweep.tex") colspec(ll cc cc) ///
    script(gen/tab_jobtype_sweep.do)
file write _t "& & \multicolumn{2}{c}{\textbf{Hierarchical}} & \multicolumn{2}{c}{\textbf{Leiden}} \\" _n
file write _t "\cmidrule(lr){3-4} \cmidrule(lr){5-6}" _n
file write _t "\textbf{Geography} & \textbf{Job type} & NMI & ARI & NMI & ARI \\" _n
file write _t "\midrule" _n
local prev ""
forvalues i = 1/`=_N' {
    local g    = geo[`i']
    local glab = cond("`g'" != "`prev'", strproper("`g'"), "")
    if "`g'" != "`prev'" & "`prev'" != "" file write _t "\addlinespace" _n
    local prev "`g'"
    local jt = jobtype[`i']
    if "`jt'" == "JT00" local jt "\textbf{JT00}\,$^\dagger$"
    else                local jt "\texttt{`jt'}"
    bt_num, value(`=NMIhc[`i']') fmt(%5.3f)
    local nh `"`r(s)'"'
    bt_num, value(`=ARIhc[`i']') fmt(%5.3f)
    local ah `"`r(s)'"'
    bt_num, value(`=NMIleiden[`i']') fmt(%5.3f)
    local nl `"`r(s)'"'
    bt_num, value(`=ARIleiden[`i']') fmt(%5.3f)
    local al `"`r(s)'"'
    bt_row, handle(_t) cells(`" "`glab'" "`jt'" "`nh'" "`ah'" "`nl'" "`al'" "')
}
bt_close, handle(_t)

display as result "  -> jobtype_sweep.tex + jobtype_sweep.csv"
