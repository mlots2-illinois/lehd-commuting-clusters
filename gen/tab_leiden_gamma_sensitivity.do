clear all

* tab_leiden_gamma_sensitivity.do -- write the Leiden resolution perturbation table (NMI and ARI against the baseline) from assess_04_sensitivity.do.
* Called by: _master.do (generate phase, do_generate switch).
* Reads:  includes _inc/booktabs.do; $tables/sensitivity_summary.dta
* Writes: $tables/leiden_gamma_sensitivity.tex; $tables/leiden_gamma_sensitivity.csv
* Notes:  Only the Leiden rows are shown, since gamma does not enter the HC run.

include "$program/_inc/booktabs.do"

local src "$tables/sensitivity_summary.dta"
capture confirm file "`src'"
if _rc {
    display as text "tab_leiden_gamma_sensitivity: missing `src' (run robust_06 + assess_04); skipping."
    exit 0
}

use "`src'", clear
keep if param == "leiden_resolution" & method == "leiden"
keep geo_level value NMI ARI

gen byte _go = cond(geo_level == "tract", 1, 2)
gsort _go value

export delimited geo_level value NMI ARI using "$tables/leiden_gamma_sensitivity.csv", replace

bt_open, handle(_t) path("$tables/leiden_gamma_sensitivity.tex") colspec(llr rr) ///
    script(gen/tab_leiden_gamma_sensitivity.do)
file write _t "Geography & Parameter & Value & NMI & ARI \\" _n
file write _t "\midrule" _n
forvalues i = 1/`=_N' {
    local g = strproper(geo_level[`i'])
    bt_num, value(`=value[`i']') fmt(%5.3f)
    local vv `"`r(s)'"'
    bt_num, value(`=NMI[`i']') fmt(%5.3f)
    local nn `"`r(s)'"'
    bt_num, value(`=ARI[`i']') fmt(%5.3f)
    local aa `"`r(s)'"'
    bt_row, handle(_t) cells(`" "`g'" "\(\gamma\)" "`vv'" "`nn'" "`aa'" "')
}
bt_close, handle(_t)

display as result "  -> leiden_gamma_sensitivity.tex + leiden_gamma_sensitivity.csv"
