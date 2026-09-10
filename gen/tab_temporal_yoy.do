clear all

* tab_temporal_yoy.do -- write the year-over-year ARI table for consecutive single-year partitions from assess_03_temporal.do as a booktabs fragment.
* Called by: _master.do (generate phase, do_generate switch).
* Reads:  includes _inc/booktabs.do; $tables/temporal_yoy.dta
* Writes: $tables/temporal_yoy.tex; $tables/temporal_yoy.csv
* Notes:  Only adjacent year pairs are kept, and only ARI is shown.

include "$program/_inc/booktabs.do"

local src "$tables/temporal_yoy.dta"
capture confirm file "`src'"
if _rc {
    display as text "tab_temporal_yoy: missing `src' (run robust_02 + assess_03); skipping."
    exit 0
}

use "`src'", clear
keep if adjacent == 1
export delimited using "$tables/temporal_yoy.csv", replace

quietly levelsof year_a if adjacent==1, local(yas)

capture program drop _cell
program define _cell, rclass
    syntax , g(string) m(string) ya(string) metric(string)
    quietly summarize `metric' if geo_level=="`g'" & method=="`m'" ///
        & year_a==`ya' & adjacent==1, meanonly
    if r(N) == 0 return local v "---"
    else         return local v = trim("`: display %5.3f r(mean)'")
end

bt_open, handle(_t) path("$tables/temporal_yoy.tex") colspec(l cc cc) ///
    script(gen/tab_temporal_yoy.do)
file write _t "& \multicolumn{2}{c}{\textbf{Tract}} & \multicolumn{2}{c}{\textbf{County}} \\" _n
file write _t "\cmidrule(lr){2-3} \cmidrule(lr){4-5}" _n
file write _t "\textbf{Transition} & HC & Leiden & HC & Leiden \\" _n
file write _t "\midrule" _n

foreach ya of local yas {
    local yb = `ya' + 1
    local cells `" "`ya'--`yb'" "'
    foreach g in tract county {
        foreach m in hc leiden {
            _cell, g("`g'") m("`m'") ya("`ya'") metric("ARI")
            local cells `"`cells' "`r(v)'""'
        }
    }
    bt_row, handle(_t) cells(`"`cells'"')
}
bt_close, handle(_t)

display as result "  -> temporal_yoy.tex + temporal_yoy.csv"
