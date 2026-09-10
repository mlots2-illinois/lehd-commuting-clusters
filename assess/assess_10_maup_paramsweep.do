clear all
cap log close
local _logdate = subinstr("$S_DATE", " ", "", .)
log using "$log/assess_10_maup_paramsweep_${run_tag}_`_logdate'.log", replace text

* assess_10_maup_paramsweep.do -- re-estimate the earnings-on-BA slope of assess_09_maup.do on every tract partition from the parameter perturbations, to show that the MAUP result does not hinge on one calibration.
* Called by: _master.do (assess phase, do_assess switch).
* Reads:  $temp/acs_tract_2023.csv (fetched by shell "$py" _py/fetch_acs.py when absent); $temp/final_assignments_<tag>_tract.dta for the baseline and the sens_<p>_m<mult> tags; $tables/maup_baseline.dta
* Writes: $tables/maup_paramsweep_<run_tag>.dta; $tables/maup_paramsweep_<run_tag>.csv; $tables/maup_paramsweep.tex (baseline run only)
* Notes:  The tag lists and parameter values are written out by hand from the calibrated tract values (alpha 0.8, cap 100 km, k 9000, gamma 0.171) rather than rebuilt from _inc/sens_grid.do. The alpha row stops at 1.25 because 1.5 x 0.8 also clips to 1 and _inc/sens_grid.do drops the duplicate.

local acs_csv "$temp/acs_tract_2023.csv"
capture confirm file "`acs_csv'"
if _rc {
    display as text "assess_10: fetching ACS tract data ..."
    shell "$py" "$program/_py/fetch_acs.py" "`acs_csv'"
}
capture confirm file "`acs_csv'"
if _rc {
    display as error "assess_10: ACS data unavailable; skipping."
    cap log close
    exit 0
}

import delimited "`acs_csv'", clear stringcols(1) varnames(1)
foreach v in agg_earn n_earn ed_tot ed_ba pop {
    capture confirm numeric variable `v'
    if _rc destring `v', replace force
}
keep if !mi(agg_earn) & !mi(n_earn) & !mi(ed_tot) & !mi(ed_ba) & !mi(pop) ///
    & n_earn > 0 & ed_tot > 0 & agg_earn > 0
keep geoid agg_earn n_earn ed_tot ed_ba pop
tempfile acsT
save "`acsT'"

tempname P
tempfile res
postfile `P' str8 param double value double mult str8 method ///
    long n_clusters double b double se double r2 using "`res'", replace

capture program drop _maup_fit
program define _maup_fit, rclass
    syntax , groupvar(string) acs(string)
    merge 1:1 geoid using "`acs'", keep(match) nogenerate
    collapse (sum) agg_earn n_earn ed_tot ed_ba, by(`groupvar')
    keep if n_earn > 0 & ed_tot > 0 & agg_earn > 0
    gen double y = agg_earn / n_earn
    gen double x = 100 * ed_ba / ed_tot
    quietly regress y x
    return scalar b  = _b[x]
    return scalar se = _se[x]
    return scalar n  = e(N)
    return scalar r2 = e(r2)
end

foreach param in alpha cap k gamma {
    if "`param'" == "alpha" {
        local tags "sens_a_m050 sens_a_m075 baseline sens_a_m125"
        local vals "0.40 0.60 0.80 1.00"
        local mlts "0.50 0.75 1.00 1.25"
        local meths "hc leiden"
    }
    if "`param'" == "cap" {
        local tags "sens_cap_m050 sens_cap_m075 baseline sens_cap_m125 sens_cap_m150"
        local vals "50 75 100 125 150"
        local mlts "0.50 0.75 1.00 1.25 1.50"
        local meths "hc leiden"
    }
    if "`param'" == "k" {
        local tags "sens_k_m050 sens_k_m075 baseline sens_k_m125 sens_k_m150"
        local vals "4500 6750 9000 11250 13500"
        local mlts "0.50 0.75 1.00 1.25 1.50"
        local meths "hc"
    }
    if "`param'" == "gamma" {
        local tags "sens_res_m050 sens_res_m075 baseline sens_res_m125 sens_res_m150"
        local vals "0.0855 0.12825 0.171 0.21375 0.2565"
        local mlts "0.50 0.75 1.00 1.25 1.50"
        local meths "leiden"
    }

    local nt : word count `tags'
    forvalues i = 1/`nt' {
        local tag : word `i' of `tags'
        local val : word `i' of `vals'
        local mlt : word `i' of `mlts'
        local af "$temp/final_assignments_`tag'_tract.dta"
        capture confirm file "`af'"
        if _rc {
            display as error "  [`param' x`mlt'] missing `tag'; skipping point."
            continue
        }
        foreach m of local meths {
            use geoid `m'_cluster using "`af'", clear
            quietly count if !mi(`m'_cluster)
            if r(N) == 0 {
                display as error "  [`param' `m' x`mlt'] no `m'_cluster; skipping."
                continue
            }
            quietly _maup_fit, groupvar(`m'_cluster) acs("`acsT'")
            post `P' ("`param'") (`val') (`mlt') ("`m'") ///
                (r(n)) (r(b)) (r(se)) (r(r2))
            display as text "  `param'=`val' (x`mlt') `m': k=" %6.0f r(n) ///
                "  b=" %7.1f r(b) " (se " %5.1f r(se) ")"
        }
    }
}
postclose `P'

use "`res'", clear
gen lo = b - invnormal(0.975) * se
gen hi = b + invnormal(0.975) * se

local b_tract .
local b_county .
capture confirm file "$tables/maup_baseline.dta"
if !_rc {
    preserve
    use "$tables/maup_baseline.dta", clear
    quietly summarize b if key == "tract",  meanonly
    local b_tract = r(mean)
    quietly summarize b if key == "county", meanonly
    local b_county = r(mean)
    restore
}
gen double ref_tract  = `b_tract'
gen double ref_county = `b_county'

label variable param      "Calibration parameter"
label variable value      "Parameter value"
label variable mult       "Multiple of chosen value"
label variable method     "Algorithm"
label variable n_clusters "Cluster count at this point"
label variable b          "OLS slope (\$/+1pp BA)"
label variable ref_tract  "Tract (finest) reference slope"
label variable ref_county "County (administrative) reference slope"
format b se lo hi ref_tract ref_county %9.1f
format r2 %5.3f

save           "$tables/maup_paramsweep_${run_tag}.dta", replace
export delimited using "$tables/maup_paramsweep_${run_tag}.csv", replace
display as result "  -> $tables/maup_paramsweep_${run_tag}.{dta,csv}"

display as text _n "{hline 74}"
display as text "MAUP robustness: slope range across each parameter sweep"
display as text "  reference: tract=" %6.1f `b_tract' "  county=" %6.1f `b_county'
display as text "{hline 74}"
preserve
collapse (min) bmin=b (max) bmax=b (min) kmin=n_clusters (max) kmax=n_clusters, ///
    by(param method)
gen byte above_county = bmin > `b_county'
list param method bmin bmax kmin kmax above_county, noobs
restore

if "${run_tag}" == "baseline" {
    capture file close _tex
    file open _tex using "$tables/maup_paramsweep.tex", write replace
    file write _tex "% Auto-generated by assess_10_maup_paramsweep.do; do not edit by hand." _n
    file write _tex "\begin{tabular}{llrrrc}" _n
    file write _tex "\toprule" _n
    file write _tex "Parameter & Method & \multicolumn{2}{c}{Slope range} & Cluster-count range & Above county? \\" _n
    file write _tex "\midrule" _n
    preserve
    collapse (min) bmin=b (max) bmax=b (min) kmin=n_clusters (max) kmax=n_clusters, ///
        by(param method)
    gen byte above_county = bmin > `b_county'
    local _N = _N
    forvalues r = 1/`_N' {
        local pp = param[`r']
        local mm = method[`r']
        local lo1 : display %9.0fc bmin[`r']
        local hi1 : display %9.0fc bmax[`r']
        local lk  : display %9.0fc kmin[`r']
        local hk  : display %9.0fc kmax[`r']
        local lo1 = subinstr(trim("`lo1'"), ",", "{,}", .)
        local hi1 = subinstr(trim("`hi1'"), ",", "{,}", .)
        local lk  = subinstr(trim("`lk'"),  ",", "{,}", .)
        local hk  = subinstr(trim("`hk'"),  ",", "{,}", .)
        local yn = cond(above_county[`r'], "yes", "no")
        file write _tex "`pp' & `mm' & `lo1' & `hi1' & `lk'--`hk' & `yn' \\" _n
    }
    restore
    file write _tex "\bottomrule" _n
    file write _tex "\end{tabular}" _n
    file close _tex
    display as result "  -> maup_paramsweep.tex"
}
else display as text "tex skipped: run_tag=${run_tag} != baseline"

display as result _n "assess_10 done."
cap log close
