clear all

* tab_dbar_sweep.do -- write the distance-cap (dbar) perturbation table from assess_04_sensitivity.do, with the chosen cap as a bold reference row.
* Called by: _master.do (generate phase, do_generate switch).
* Reads:  includes _inc/booktabs.do and _inc/get_chosen.do; $tables/sensitivity_summary.dta; $temp/chosen_params_<lvl>.dta
* Writes: $tables/dbar_sweep.tex; $tables/dbar_sweep.csv
* Notes:  The baseline row is not in sensitivity_summary.dta, so we add it with NMI = ARI = 1 when no perturbation lands on the chosen value.

include "$program/_inc/booktabs.do"
include "$program/_inc/get_chosen.do"

local src "$tables/sensitivity_summary.dta"
capture confirm file "`src'"
if _rc {
    display as text "tab_dbar_sweep: missing `src' (run robust_06 + assess_04); skipping."
    exit 0
}

local c_t = 100
local c_c = 1500
get_chosen, lvl(tract) param(cap_km)
if r(value) < . local c_t = r(value)
get_chosen, lvl(county) param(cap_km)
if r(value) < . local c_c = r(value)

use "`src'", clear
keep if param == "cap_km"
keep geo_level value method NMI ARI
reshape wide NMI ARI, i(geo_level value) j(method) string

gen byte _isbase = 0

foreach lvl of global levels {
    local bval = cond("`lvl'" == "tract", `c_t', `c_c')
    quietly count if geo_level == "`lvl'" & abs(value - `bval') < 1e-6*max(1, abs(`bval'))
    if r(N) == 0 {
        local n0 = _N + 1
        quietly set obs `n0'
        quietly replace geo_level = "`lvl'" in `n0'
        quietly replace value     = `bval'  in `n0'
        foreach v in NMIhc ARIhc NMIleiden ARIleiden {
            quietly replace `v' = 1 in `n0'
        }
        quietly replace _isbase = 1 in `n0'
    }
}
replace _isbase = 1 if (geo_level=="tract"  & abs(value-`c_t') < 1e-6*max(1, abs(`c_t'))) ///
    | (geo_level=="county" & abs(value-`c_c') < 1e-6*max(1, abs(`c_c')))

gen byte _go = cond(geo_level == "tract", 1, 2)
gsort _go value

export delimited geo_level value NMIhc ARIhc NMIleiden ARIleiden ///
    using "$tables/dbar_sweep.csv", replace

bt_open, handle(_t) path("$tables/dbar_sweep.tex") colspec(ll cc cc) ///
    script(gen/tab_dbar_sweep.do)
file write _t "& & \multicolumn{2}{c}{\textbf{Hierarchical}} & \multicolumn{2}{c}{\textbf{Leiden}} \\" _n
file write _t "\cmidrule(lr){3-4} \cmidrule(lr){5-6}" _n
file write _t "\textbf{Geography} & \(\bar{d}\) (km) & NMI & ARI & NMI & ARI \\" _n
file write _t "\midrule" _n

local prev ""
forvalues i = 1/`=_N' {
    local g    = geo_level[`i']
    local glab = cond("`g'" != "`prev'", strproper("`g'"), "")
    if "`g'" != "`prev'" & "`prev'" != "" file write _t "\addlinespace" _n
    local prev "`g'"

    if _isbase[`i'] bt_num, value(`=value[`i']') fmt(%5.0f) bold dagger
    else            bt_num, value(`=value[`i']') fmt(%5.0f)
    local vv `"`r(s)'"'
    bt_num, value(`=NMIhc[`i']') fmt(%5.3f)
    local nh `"`r(s)'"'
    bt_num, value(`=ARIhc[`i']') fmt(%5.3f)
    local ah `"`r(s)'"'
    bt_num, value(`=NMIleiden[`i']') fmt(%5.3f)
    local nl `"`r(s)'"'
    bt_num, value(`=ARIleiden[`i']') fmt(%5.3f)
    local al `"`r(s)'"'
    bt_row, handle(_t) cells(`" "`glab'" "`vv'" "`nh'" "`ah'" "`nl'" "`al'" "')
}
bt_close, handle(_t)

display as result "  -> dbar_sweep.tex + dbar_sweep.csv"
