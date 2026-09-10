clear all
cap log close
local _logdate = subinstr("$S_DATE", " ", "", .)
log using "$log/report_01_descriptive_table_${run_tag}_`_logdate'.log", replace text

* report_01_descriptive_table.do -- summarize cluster count, units per cluster and population per cluster for HC and Leiden at each geography, and write the descriptive table.
* Called by: _master.do (generate phase, do_generate switch).
* Reads:  includes _inc/booktabs.do; $temp/final_assignments_<run_tag>_<lvl>.dta; $temp/centroids_tract.dta
* Writes: $tables/descriptive_<method>_<lvl>.dta and .csv; $tables/descriptive_clusters.tex
* Notes:  County population is the sum of tract populations from centroids_tract.dta, so both geographies rest on the same 2020 counts. Units without a population record count as zero rather than dropping.

include "$program/_inc/booktabs.do"

capture program drop _build_descriptive
program define _build_descriptive
    syntax , src(string) clusvar(name) outfile(string) pops(string) method(string)

    use "`src'", clear
    quietly count
    local _n_assigned = r(N)
    merge m:1 geoid using "`pops'", keep(master match) keepusing(pop) nogenerate
    quietly count if mi(pop)
    local _n_nopop = r(N)
    if `_n_nopop' > 0 {
        display as text "  [`method'] " `_n_nopop' " of " `_n_assigned' ///
            " assigned units had no pop record -> pop=0"
    }
    replace pop = 0 if mi(pop)

    egen long _global_cluster = group(`clusvar')

    bysort _global_cluster: egen long   _ntr = count(geoid)
    bysort _global_cluster: egen double _pop = total(pop)
    bysort _global_cluster: keep if _n == 1

    gcollapse                               ///
        (count) n_clusters = _ntr           ///
        (mean)  mean_pop   = _pop           ///
        (sd)    sd_pop     = _pop           ///
        (p50)   med_pop    = _pop           ///
        (min)   min_pop    = _pop           ///
        (max)   max_pop    = _pop           ///
        (mean)  mean_units = _ntr           ///
        (sd)    sd_units   = _ntr           ///
        (p50)   med_units  = _ntr           ///
        (min)   min_units  = _ntr           ///
        (max)   max_units  = _ntr

    gen str8 method = "`method'"
    order method
    format mean_pop sd_pop med_pop min_pop max_pop %10.0fc
    format mean_units sd_units med_units min_units max_units %6.2f

    save "`outfile'.dta", replace
    export delimited using "`outfile'.csv", replace
end

foreach lvl of global levels {

    display as text _n "{hline 70}"
    display as text "===== report_01 LEVEL: `lvl' ====="
    display as text "{hline 70}"

    capture confirm file "$temp/final_assignments_${run_tag}_`lvl'.dta"
    if _rc {
        display as text "report_01: no final_assignments for `lvl'; skipping level."
        continue
    }
    capture confirm file "$temp/centroids_tract.dta"
    if _rc {
        display as text "report_01: $temp/centroids_tract.dta not found (needed for the pop lookup); skipping level `lvl'."
        continue
    }

    if "`lvl'" == "tract" {
        use geoid pop using "$temp/centroids_tract.dta", clear
    }
    else {
        use geoid pop using "$temp/centroids_tract.dta", clear
        gen str5 county_geoid = substr(geoid, 1, 5)
        gcollapse (sum) pop, by(county_geoid)
        rename county_geoid geoid
    }
    tempfile pops
    save `pops'

    local src "$temp/final_assignments_${run_tag}_`lvl'.dta"

    _build_descriptive,                                ///
        src("`src'") clusvar(hc_cluster)               ///
        outfile("$tables/descriptive_hc_`lvl'")        ///
        pops("`pops'") method("HC")

    display as text _n "HC descriptive (`lvl'):"
    list, noobs

    _build_descriptive,                                ///
        src("`src'") clusvar(leiden_cluster)           ///
        outfile("$tables/descriptive_leiden_`lvl'")    ///
        pops("`pops'") method("Leiden")

    display as text _n "Leiden descriptive (`lvl'):"
    list, noobs

    display as result "report_01 (`lvl') done."
}

tempfile _desc_acc
local _have 0
foreach lvl in tract county {
    foreach m in hc leiden {
        capture confirm file "$tables/descriptive_`m'_`lvl'.dta"
        if _rc continue
        use "$tables/descriptive_`m'_`lvl'.dta", clear
        gen str6 geo_level = "`lvl'"
        gen byte _ordg = cond("`lvl'" == "tract", 1, 2)
        gen byte _ordm = cond("`m'"   == "hc",    1, 2)
        if `_have' append using "`_desc_acc'"
        save "`_desc_acc'", replace
        local _have 1
    }
}
if !`_have' {
    display as text "report_01: no per-level descriptive files exist; descriptive_clusters.tex not written."
    cap log close
    exit 0
}
use "`_desc_acc'", clear
sort _ordg _ordm

bt_open, handle(_tex) path("$tables/descriptive_clusters.tex") colspec(llrrrrrrr) ///
    script(report/report_01_descriptive_table.do)
file write _tex "Geography & Method & \(n_{\text{clusters}}\) & Mean units & SD units & Mean pop. & SD pop. & Median pop. & Max pop. \\" _n
file write _tex "\midrule" _n
local _N = _N
forvalues r = 1/`_N' {
    local geo = strproper(geo_level[`r'])
    local mth = method[`r']
    bt_num, value(`=n_clusters[`r']') fmt(%9.0fc)
    local nc `"`r(s)'"'
    bt_num, value(`=mean_units[`r']') fmt(%6.2f)
    local mu `"`r(s)'"'
    bt_num, value(`=sd_units[`r']') fmt(%6.2f)
    local su `"`r(s)'"'
    bt_num, value(`=mean_pop[`r']') fmt(%12.0fc)
    local mp `"`r(s)'"'
    bt_num, value(`=sd_pop[`r']') fmt(%12.0fc)
    local sp `"`r(s)'"'
    bt_num, value(`=med_pop[`r']') fmt(%12.0fc)
    local dp `"`r(s)'"'
    bt_num, value(`=max_pop[`r']') fmt(%12.0fc)
    local xp `"`r(s)'"'
    bt_row, handle(_tex) cells(`" "`geo'" "`mth'" "`nc'" "`mu'" "`su'" "`mp'" "`sp'" "`dp'" "`xp'" "')
}
bt_close, handle(_tex)
display as result "  -> descriptive_clusters.tex (LaTeX fragment for tab:descriptive_clusters)"

display as result _n "report_01 done."

cap log close
