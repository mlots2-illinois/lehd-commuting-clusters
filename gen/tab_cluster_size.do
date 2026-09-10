clear all

* tab_cluster_size.do -- write the cluster-size summary table (count, size quantiles, maximum, share of singletons) by geography and method.
* Called by: _master.do (generate phase, do_generate switch).
* Reads:  includes _inc/booktabs.do and _inc/cluster_sizes.do; $temp/final_assignments_<run_tag>_<lvl>.dta
* Writes: $tables/cluster_size_summary.tex; $tables/cluster_size_summary.csv

include "$program/_inc/booktabs.do"
include "$program/_inc/cluster_sizes.do"

tempfile acc
local first 1

foreach lvl of global levels {
    capture confirm file "$temp/final_assignments_${run_tag}_`lvl'.dta"
    if _rc {
        display as text "  no final_assignments for `lvl'; skipping."
        continue
    }
    use geoid hc_cluster leiden_cluster using ///
        "$temp/final_assignments_${run_tag}_`lvl'.dta", clear

    cluster_sizes

    gen byte _singleton = (size == 1)
    gen str6 geo_level = "`lvl'"
    gcollapse (count) n_clusters = size                          ///
        (p10)  p10  = size (p25) p25 = size (p50) p50 = size ///
        (p75)  p75  = size (p90) p90 = size (p99) p99 = size ///
        (max)  max_size = size                             ///
        (mean) pct_singleton = _singleton,                 ///
        by(geo_level method)
    replace pct_singleton = 100 * pct_singleton

    if `first' {
        save `acc', replace
        local first 0
    }
    else {
        append using `acc'
        save `acc', replace
    }
}

capture confirm file "`acc'"
if _rc {
    display as text "tab_cluster_size: nothing to write; skipping."
    exit 0
}

use `acc', clear
gen byte _ord = cond(geo_level == "tract", 1, 2)
gen byte _ordm = cond(method == "HC", 1, 2)
sort _ord _ordm
drop _ord _ordm
order geo_level method n_clusters p10 p25 p50 p75 p90 p99 max_size pct_singleton
export delimited using "$tables/cluster_size_summary.csv", replace

bt_open, handle(_t) path("$tables/cluster_size_summary.tex") colspec(lllrrrrrrrr) ///
    script(gen/tab_cluster_size.do)
file write _t "          &        & \(n_{\text{clust}}\) & \multicolumn{6}{c}{Cluster size quantiles (units per cluster)} & Max & \% \\" _n
file write _t "Geography & Method &                    & p10 & p25 & p50 & p75 & p90 & p99 & size & singleton \\" _n
file write _t "\midrule" _n
forvalues r = 1/`=_N' {
    local geo = strproper(geo_level[`r'])
    local mth = method[`r']
    bt_num, value(`=n_clusters[`r']') fmt(%9.0fc)
    local nc `"`r(s)'"'
    foreach q in p10 p25 p50 p75 p90 p99 max_size {
        bt_num, value(`=`q'[`r']') fmt(%9.0fc)
        local _`q' `"`r(s)'"'
    }
    bt_num, value(`=pct_singleton[`r']') fmt(%4.1f)
    local ps `"`r(s)'"'
    bt_row, handle(_t) cells(`" "`geo'" "`mth'" "`nc'" "`_p10'" "`_p25'" "`_p50'" "`_p75'" "`_p90'" "`_p99'" "`_max_size'" "`ps'" "')
}
bt_close, handle(_t)
display as result "  -> cluster_size_summary.tex + cluster_size_summary.csv"
