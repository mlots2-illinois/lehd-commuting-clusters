clear all
cap log close
local _logdate = subinstr("$S_DATE", " ", "", .)
log using "$log/assess_04_sensitivity_${run_tag}_`_logdate'.log", replace text

* assess_04_sensitivity.do -- compare every parameter-perturbation rerun from robust_06_params.do with the baseline partition and summarize the agreement.
* Called by: _master.do (assess phase, do_assess switch).
* Reads:  includes _inc/compare_partitions.do and _inc/sens_grid.do; $temp/final_assignments_baseline_<lvl>.dta; $temp/final_assignments_sens_<p>_m<mult>_<lvl>.dta
* Writes: $tables/sensitivity_summary.dta; $tables/sensitivity_summary.tex
* Notes:  The list of perturbation tags and values is rebuilt from _inc/sens_grid.do rather than read from disk, so the same grid drives the reruns and this summary. Runs that are missing on disk are listed and skipped; the alpha endpoint runs from robust_08_alpha_full.do are not in the grid and are picked up by tab_alpha_sweep.do instead.

include "$program/_inc/compare_partitions.do"

include "$program/_inc/sens_grid.do"

local methods "hc leiden"

tempfile out
local first 1

foreach lvl of global levels {

    display as text _n "{hline 70}"
    display as text "===== assess_04 LEVEL: `lvl' ====="
    display as text "{hline 70}"

    capture confirm file "$temp/final_assignments_baseline_`lvl'.dta"
    if _rc {
        display as error "Missing baseline run for `lvl'; skipping."
        continue
    }

    local n_specs : word count `sens_tags'
    local present_idx ""
    local missing ""
    forvalues i = 1/`n_specs' {
        local tag : word `i' of `sens_tags'
        capture confirm file "$temp/final_assignments_`tag'_`lvl'.dta"
        if _rc {
            local missing "`missing' `tag'"
        }
        else {
            local present_idx "`present_idx' `i'"
        }
    }
    if "`missing'" != "" {
        display as text "  Skipping (no file at `lvl'):`missing'"
    }
    if "`present_idx'" == "" {
        display as text "  No perturbation runs at `lvl'; skipping."
        continue
    }

    use geoid hc_cluster leiden_cluster ///
        using "$temp/final_assignments_baseline_`lvl'.dta", clear
    rename hc_cluster     hc_base
    rename leiden_cluster le_base
    quietly count
    local n_base = r(N)
    tempfile base
    save `base'

    foreach i of local present_idx {
        local tag    : word `i' of `sens_tags'
        local pname  : word `i' of `sens_params'
        local v_t    : word `i' of `sens_v_t'
        local v_c    : word `i' of `sens_v_c'
        local pval   = cond("`lvl'" == "tract", "`v_t'", "`v_c'")

        use geoid hc_cluster leiden_cluster ///
            using "$temp/final_assignments_`tag'_`lvl'.dta", clear
        rename hc_cluster     hc_pert
        rename leiden_cluster le_pert
        quietly count
        local n_pert = r(N)
        merge 1:1 geoid using `base', keep(match) nogenerate
        quietly count
        local n_both      = r(N)
        local n_dropped   = `n_pert' - `n_both'
        local n_dropped_b = `n_base' - `n_both'
        if `n_dropped' > 0 {
            display as text "  [`lvl' `tag'] " `n_dropped' " of " `n_pert' ///
                " perturbation-units not in baseline (compared on " ///
                `n_both' "-unit intersection)"
        }
        if `n_dropped_b' > 0 {
            display as text "  [`lvl' `tag'] " `n_dropped_b' " of " `n_base' ///
                " baseline-units not in perturbation (compared on " ///
                `n_both' "-unit intersection)"
        }

        foreach m of local methods {
            local va = cond("`m'" == "hc", "hc_base", "le_base")
            local vb = cond("`m'" == "hc", "hc_pert", "le_pert")

            compare_partitions, a(`va') b(`vb')

            preserve
            clear
            quietly set obs 1
            gen str6   geo_level = "`lvl'"
            gen str20  run_tag   = "`tag'"
            gen str20  param     = "`pname'"
            gen double value     = `pval'
            gen str8   method    = "`m'"
            gen long   n         = `r(N)'
            gen long   n_dropped = `n_dropped'
            gen double NMI       = `r(NMI)'
            gen double ARI       = `r(ARI)'
            gen double AMI       = `r(AMI)'

            if `first' {
                save `out', replace
                local first 0
            }
            else {
                append using `out'
                save `out', replace
            }
            restore
        }

        display as text "  `tag' (`lvl'): done"
    }
}

capture confirm file "`out'"
if _rc {
    display as error "No sensitivity comparisons produced."
    display as error "  robust_06_params.do produces them (do_robust_params in _master.do)."
    cap log close
    exit 0
}

use `out', clear
sort geo_level param value method

label variable geo_level "Geography (tract or county)"
label variable run_tag   "Perturbation run tag"
label variable param     "Perturbed parameter"
label variable value     "Perturbed value"
label variable method    "Algorithm"
label variable n         "Units compared (intersection with baseline)"
label variable n_dropped "Perturbation units absent from baseline"
label variable NMI       "NMI vs. baseline (arithmetic; NOT chance-corrected)"
label variable ARI       "ARI vs. baseline (chance-corrected)"
label variable AMI       "AMI vs. baseline (chance-corrected)"

format NMI ARI AMI %7.4f
save "$tables/sensitivity_summary.dta", replace

display as text _n "{hline 70}"
display as text "Sensitivity vs. baseline"
display as text "{hline 70}"
list geo_level param value method n AMI ARI, noobs sepby(geo_level param)

capture file close _tex
file open _tex using "$tables/sensitivity_summary.tex", write replace
file write _tex "% Auto-generated by assess_04_sensitivity.do; do not edit by hand." _n
file write _tex "\begin{tabular}{llrr}" _n
file write _tex "\toprule" _n
file write _tex "Geography & Perturbation & ARI (HC) & ARI (Leiden) \\" _n
file write _tex "\midrule" _n
foreach g of global levels {
    local glab = strproper("`g'")
    local wrote_g 0
    foreach tag of local sens_tags {
        quietly count if geo_level == "`g'" & run_tag == "`tag'"
        if r(N) == 0 continue
        quietly levelsof param if geo_level=="`g'" & run_tag=="`tag'", local(pn) clean
        quietly summarize value if geo_level=="`g'" & run_tag=="`tag'", meanonly
        local pv = r(mean)
        if "`pn'" == "alpha"             local plab "\(\alpha\) = `: display %4.2f `pv''"
        else if "`pn'" == "cap_km"       local plab "\(\bar{d}\) = `: display %5.0f `pv'' km"
        else if "`pn'" == "hc_k_national" local plab "\(k\) = `: display %9.0fc `pv''"
        else if "`pn'" == "leiden_resolution" local plab "\(\gamma\) = `: display %5.2f `pv''"
        else local plab "`pn' = `pv'"
        local plab = subinstr("`plab'", ",", "{,}", .)
        quietly summarize ARI if geo_level=="`g'" & run_tag=="`tag'" & method=="hc", meanonly
        local ah = cond(r(N)==0, "---", "`: display %5.3f r(mean)'")
        quietly summarize ARI if geo_level=="`g'" & run_tag=="`tag'" & method=="leiden", meanonly
        local al = cond(r(N)==0, "---", "`: display %5.3f r(mean)'")
        local gcell = cond(`wrote_g'==0, "`glab'", "")
        file write _tex "`gcell' & `plab' & `=trim("`ah'")' & `=trim("`al'")' \\" _n
        local wrote_g 1
    }
}
file write _tex "\bottomrule" _n
file write _tex "\end{tabular}" _n
file close _tex
display as result "  -> sensitivity_summary.tex (LaTeX fragment for tab:sensitivity_summary)"

display as result _n "assess_04 done."

cap log close
