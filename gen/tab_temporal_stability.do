clear all

* tab_temporal_stability.do -- write the window-stability table (5-year versus 3-year versus 1-year) from assess_03_temporal.do as a booktabs fragment.
* Called by: _master.do (generate phase, do_generate switch).
* Reads:  includes _inc/booktabs.do; $tables/temporal_agreement.dta
* Writes: $tables/temporal_stability.tex; $tables/temporal_stability.csv
* Notes:  Reads the agreement statistics that assess_03_temporal.do saved and lays them out as a booktabs fragment for the manuscript.

include "$program/_inc/booktabs.do"

local src "$tables/temporal_agreement.dta"
capture confirm file "`src'"
if _rc {
    display as text "tab_temporal_stability: missing `src' (run robust_01/02 + assess_03); skipping."
    exit 0
}

use "`src'", clear
export delimited using "$tables/temporal_stability.csv", replace

capture program drop _cell
program define _cell, rclass
    syntax , g(string) m(string) pa(string) pb(string) metric(string)
    quietly summarize `metric' if geo_level=="`g'" & method=="`m'" ///
        & period_a=="`pa'" & period_b=="`pb'", meanonly
    if r(N) == 0 return local v "---"
    else         return local v = trim("`: display %5.3f r(mean)'")
end

bt_open, handle(_t) path("$tables/temporal_stability.tex") colspec(ll ccc ccc) ///
    script(gen/tab_temporal_stability.do)
file write _t "& & \multicolumn{3}{c}{\textbf{NMI}} & \multicolumn{3}{c}{\textbf{ARI}} \\" _n
file write _t "\cmidrule(lr){3-5} \cmidrule(lr){6-8}" _n
file write _t "\textbf{Geography} & \textbf{Method} & 5yr v 3yr & 5yr v 1yr & 3yr v 1yr & 5yr v 3yr & 5yr v 1yr & 3yr v 1yr \\" _n
file write _t "\midrule" _n

local first 1
foreach g of global levels {
    local glab = strproper("`g'")
    if !`first' file write _t "\addlinespace" _n
    local first 0
    local wrote_g 0
    foreach m in hc leiden {
        local mlab = cond("`m'"=="hc", "HC", "Leiden")
        local gcell = cond(`wrote_g'==0, "`glab'", "")
        local cells `" "`gcell'" "`mlab'" "'
        foreach metric in NMI ARI {
            foreach pair in "5yr 3yr" "5yr 1yr" "3yr 1yr" {
                local pa : word 1 of `pair'
                local pb : word 2 of `pair'
                _cell, g("`g'") m("`m'") pa("`pa'") pb("`pb'") metric("`metric'")
                local cells `"`cells' "`r(v)'""'
            }
        }
        bt_row, handle(_t) cells(`"`cells'"')
        local wrote_g 1
    }
}
bt_close, handle(_t)

display as result "  -> temporal_stability.tex + temporal_stability.csv"
