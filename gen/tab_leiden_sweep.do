clear all

* tab_leiden_sweep.do -- write the Leiden resolution sweep table from calib_02_resolution.do, thinned to the display step, with the chosen gamma in bold.
* Called by: _master.do (generate phase, do_generate switch).
* Reads:  includes _inc/booktabs.do and _inc/get_chosen.do; $tables/leiden_resolution_sweep_<lvl>.csv; $temp/chosen_params_<lvl>.dta
* Writes: $tables/leiden_sweep.tex; $tables/leiden_sweep.csv
* Notes:  Rows are kept where the resolution is a multiple of $leiden_sweep_step, plus the row nearest the chosen gamma; a dagger marks that row when the chosen value was not itself swept. The CSV keeps every row.

include "$program/_inc/booktabs.do"
include "$program/_inc/get_chosen.do"

tempfile comb
local first 1
foreach lvl of global levels {
    local csv "$tables/leiden_resolution_sweep_`lvl'.csv"
    capture confirm file "`csv'"
    if _rc continue
    preserve
    import delimited using "`csv'", varnames(1) clear
    gen str6 geo_level = "`lvl'"
    if `first' {
        save `comb', replace
        local first 0
    }
    else {
        append using `comb'
        save `comb', replace
    }
    restore
}
capture confirm file "`comb'"
if _rc {
    display as text "tab_leiden_sweep: no leiden_resolution_sweep_*.csv found (run calib_02); skipping."
    exit 0
}
use `comb', clear
order geo_level resolution n_clusters_est_national mean_stability ///
    mean_largest_frac mean_singleton_frac mean_modularity
export delimited using "$tables/leiden_sweep.csv", replace

bt_open, handle(_t) path("$tables/leiden_sweep.tex") colspec(rrrrrr) ///
    script(gen/tab_leiden_sweep.do)
file write _t "\(\gamma\) & \(\hat{n}_{\text{nat}}\) & \(\bar{s}\) & \(\bar{f}_{\max}\) & \(\bar{f}_{1}\) & \(\bar{Q}\) \\" _n

foreach lvl of global levels {
    quietly count if geo_level == "`lvl'"
    if r(N) == 0 continue

    get_chosen, lvl(`lvl') param(leiden_resolution)
    local chosen = r(value)

    local plab = cond("`lvl'"=="tract", "Tract sweep", "County sweep")
    file write _t "\midrule" _n
    file write _t "\multicolumn{6}{c}{\textit{`plab'}} \\" _n

    preserve
    keep if geo_level == "`lvl'"
    sort resolution

    local pstep = 0
    capture local pstep = $leiden_sweep_step

    gen double _dev = abs(resolution - `chosen')
    quietly summarize _dev, meanonly
    local dev_min = r(min)
    local inexact = (`chosen' < .) & (`dev_min' > 1e-6*max(1, abs(`chosen')))
    if `inexact' display as text "  tab_leiden_sweep: chosen gamma " ///
        "`chosen' (`lvl') is not a swept value; marking the nearest row."

    if `pstep' > 0 {
        gen byte _print = ///
            abs(resolution/`pstep' - round(resolution/`pstep')) < 1e-4
        replace _print = 1 if (`chosen' < .) & (_dev <= `dev_min' + 1e-9)
        keep if _print
    }

    replace mean_modularity = 0 if abs(mean_modularity) < 0.0005

    forvalues i = 1/`=_N' {
        local isch = (`chosen' < .) & (_dev[`i'] <= `dev_min' + 1e-9)
        local bopt = cond(`isch', "bold", "")
        local dopt = cond(`isch' & `inexact', "dagger", "")
        bt_num, value(`=resolution[`i']') fmt(%5.3f) `bopt' `dopt'
        local g `"`r(s)'"'
        bt_num, value(`=n_clusters_est_national[`i']') fmt(%12.0fc) `bopt'
        local nn `"`r(s)'"'
        bt_num, value(`=mean_stability[`i']') fmt(%5.3f) `bopt'
        local sb `"`r(s)'"'
        bt_num, value(`=mean_largest_frac[`i']') fmt(%5.3f) `bopt'
        local fm `"`r(s)'"'
        bt_num, value(`=mean_singleton_frac[`i']') fmt(%5.3f) `bopt'
        local f1 `"`r(s)'"'
        bt_num, value(`=mean_modularity[`i']') fmt(%5.3f) `bopt'
        local qb `"`r(s)'"'
        bt_row, handle(_t) cells(`" "`g'" "`nn'" "`sb'" "`fm'" "`f1'" "`qb'" "')
    }
    restore
}
bt_close, handle(_t)

display as result "  -> leiden_sweep.tex + leiden_sweep.csv"
