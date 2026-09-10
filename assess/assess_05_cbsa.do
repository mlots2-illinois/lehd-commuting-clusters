clear all
cap log close
local _logdate = subinstr("$S_DATE", " ", "", .)
log using "$log/assess_05_cbsa_${run_tag}_`_logdate'.log", replace text

* assess_05_cbsa.do -- measure how well each partition preserves core-based statistical areas (CBSA), using the share of a CBSA's units that fall in its largest cluster.
* Called by: _master.do (assess phase, do_assess switch).
* Reads:  includes _inc/maxshare_purity.do; $temp/final_assignments_<run_tag>_<lvl>.dta; $temp/cbsa/tract_to_cbsa.dta; $temp/cbsa/county_to_cbsa.dta
* Writes: $tables/cbsa_preservation_<run_tag>.dta; $tables/cbsa_preservation_<run_tag>.csv; $tables/cbsa_preservation.tex (baseline run only)
* Notes:  Units outside any CBSA drop at the merge, so purity describes metropolitan and micropolitan territory only. We also carry the inverse-cluster-count purity as a check on the max-share measure, but only the max-share version reaches the LaTeX table.

include "$program/_inc/maxshare_purity.do"

local xwalk_dir "$temp/cbsa"

tempfile _acc
local _have 0

foreach lvl of global levels {

    display as text _n "{hline 70}"
    display as text "===== assess_05 LEVEL: `lvl' ====="
    display as text "{hline 70}"

    local xw = cond("`lvl'" == "tract", ///
        "`xwalk_dir'/tract_to_cbsa.dta", "`xwalk_dir'/county_to_cbsa.dta")

    capture confirm file "`xw'"
    if _rc {
        display as error "Missing CBSA crosswalk for `lvl': `xw'; skipping."
        continue
    }
    capture confirm file "$temp/final_assignments_${run_tag}_`lvl'.dta"
    if _rc {
        display as error "Missing final_assignments for `lvl'; skipping."
        continue
    }

    foreach m in hc leiden {
        local mlab = cond("`m'" == "hc", "HC", "Leiden")

        use geoid `m'_cluster ///
            using "$temp/final_assignments_${run_tag}_`lvl'.dta", clear
        rename `m'_cluster clus

        quietly count
        local n_pre = r(N)
        merge m:1 geoid using "`xw'", keep(match) keepusing(cbsa) nogenerate
        quietly count
        local n_matched = r(N)
        display as text "  CBSA merge (`lvl' `mlab'): `n_matched' of `n_pre' units matched (" ///
            %5.1f 100 * `n_matched' / `n_pre' "%)"
        if `n_matched' == 0 {
            display as error "assess_05: no units matched the CBSA crosswalk;" ///
                " check the crosswalk vintage (geoids in `xw' must match the" ///
                " LODES tract/county definitions)."
        }
        assert `n_matched' > 0

        maxshare_purity, refvar(cbsa) clustervar(clus) gen(purity) altgen(purity_alt)

        quietly count
        local n_cbsa = r(N)

        gen byte _full = (purity == 1)
        gen byte _frag = (purity < 0.5)
        quietly summarize purity, detail
        local mean = r(mean)
        local med  = r(p50)
        local p10  = r(p10)
        local p25  = r(p25)
        local p75  = r(p75)
        local p90  = r(p90)
        quietly summarize _full
        local pct_full = r(mean) * 100
        quietly summarize _frag
        local pct_frag = r(mean) * 100
        quietly summarize purity_alt
        local mean_alt = r(mean)

        display as text "  `mlab' (`lvl'): n_cbsa=`n_cbsa'  mean=" %5.3f `mean' ///
            "  med=" %5.3f `med' "  %fully=" %4.1f `pct_full' ///
            "  %highly_frag=" %4.1f `pct_frag' "  [alt mean=" %5.3f `mean_alt' "]"

        clear
        set obs 1
        gen str6   geo_level    = "`lvl'"
        gen str8   method       = "`mlab'"
        gen long   n_cbsa       = `n_cbsa'
        gen double mean_purity  = `mean'
        gen double med_purity   = `med'
        gen double p10          = `p10'
        gen double p25          = `p25'
        gen double p75          = `p75'
        gen double p90          = `p90'
        gen double pct_fully    = `pct_full'
        gen double pct_highly   = `pct_frag'
        gen double mean_purity_alt = `mean_alt'

        if `_have' append using "`_acc'"
        save "`_acc'", replace
        local _have 1
    }
}

capture confirm file "`_acc'"
if _rc {
    display as error "No CBSA preservation computed (missing crosswalks or assignments)."
    cap log close
    exit 0
}

use "`_acc'", clear
gen byte _ordg = cond(geo_level == "tract", 1, 2)
gen byte _ordm = cond(method == "HC", 1, 2)
sort _ordg _ordm
drop _ordg _ordm

label variable geo_level       "Geography (tract or county)"
label variable method          "Algorithm"
label variable n_cbsa          "CBSAs with units in the assignment"
label variable mean_purity     "Mean per-CBSA max-share purity"
label variable med_purity      "Median per-CBSA purity"
label variable pct_fully       "% CBSAs wholly in one cluster (purity = 1)"
label variable pct_highly      "% CBSAs highly fragmented (purity < 0.5)"
label variable mean_purity_alt "Mean per-CBSA inverse-cluster-count purity (alt)"
format mean_purity med_purity p10 p25 p75 p90 mean_purity_alt %5.3f
format pct_fully pct_highly %5.1f

save           "$tables/cbsa_preservation_${run_tag}.dta", replace
export delimited using "$tables/cbsa_preservation_${run_tag}.csv", replace
display as result "  -> $tables/cbsa_preservation_${run_tag}.{csv,dta}"

display as text _n "{hline 70}"
display as text "CBSA preservation (max-share purity)"
display as text "{hline 70}"
list geo_level method n_cbsa mean_purity med_purity pct_fully pct_highly, noobs

if "${run_tag}" == "baseline" {
    capture file close _tex
    file open _tex using "$tables/cbsa_preservation.tex", write replace
    file write _tex "% Auto-generated by assess_05_cbsa.do; do not edit by hand." _n
    file write _tex "\begin{tabular}{llrrrrrrrr}" _n
    file write _tex "\toprule" _n
    file write _tex "          &        & \multicolumn{6}{c}{Per-CBSA cluster purity} & \% fully  & \% highly \\" _n
    file write _tex "Geography & Method & Mean & Med & p10 & p25 & p75 & p90          & preserved & fragmented \\" _n
    file write _tex "\midrule" _n
    local _N = _N
    forvalues r = 1/`_N' {
        local geo = strproper(geo_level[`r'])
        local mth = method[`r']
        local me  : display %5.3f mean_purity[`r']
        local md  : display %5.3f med_purity[`r']
        local a   : display %5.3f p10[`r']
        local b   : display %5.3f p25[`r']
        local c   : display %5.3f p75[`r']
        local d   : display %5.3f p90[`r']
        local pf  : display %4.1f pct_fully[`r']
        local pg  : display %4.1f pct_highly[`r']
        file write _tex "`geo' & `mth' & `=trim("`me'")' & `=trim("`md'")' & `=trim("`a'")' & `=trim("`b'")' & `=trim("`c'")' & `=trim("`d'")' & `=trim("`pf'")' & `=trim("`pg'")' \\" _n
    }
    file write _tex "\bottomrule" _n
    file write _tex "\end{tabular}" _n
    file close _tex
    display as result "  -> cbsa_preservation.tex (LaTeX fragment for tab:cbsa_preservation)"
}
else display as text "tex skipped: run_tag=${run_tag} != baseline"

display as result _n "assess_05 done."

cap log close
