clear all
cap log close
local _logdate = subinstr("$S_DATE", " ", "", .)
log using "$log/assess_09_maup_${run_tag}_`_logdate'.log", replace text

* assess_09_maup.do -- estimate the slope of mean earnings on percent BA at ten areal units, from nation to tract, to place the commuting clusters on the modifiable areal unit problem (MAUP) ladder.
* Called by: _master.do (assess phase, do_assess switch).
* Reads:  includes _inc/census_regions.do; $temp/acs_tract_2023.csv (fetched by shell "$py" _py/fetch_acs.py when absent); $temp/final_assignments_<run_tag>_tract.dta; $temp/final_assignments_<run_tag>_county.dta
* Writes: $tables/maup_<run_tag>.dta; $tables/maup_<run_tag>.csv; $tables/maup.tex (baseline run only)
* Notes:  Every unit on the ladder is built by summing the same complete-case tract records, so the estimates differ only by aggregation. We report unweighted and population-weighted slopes because the weighted one is what an analyst would usually run. The ACS fetch needs network access; without the CSV the script exits cleanly.

include "$program/_inc/census_regions.do"

local acs_csv "$temp/acs_tract_2023.csv"
capture confirm file "`acs_csv'"
if _rc {
    display as text "assess_09: fetching ACS tract data ..."
    shell "$py" "$program/_py/fetch_acs.py" "`acs_csv'"
}
capture confirm file "`acs_csv'"
if _rc {
    display as error "assess_09: `acs_csv' unavailable (ACS fetch failed); skipping."
    cap log close
    exit 0
}

capture confirm file "$temp/final_assignments_${run_tag}_tract.dta"
if _rc {
    display as error "assess_09: tract assignments missing; run the tract pipeline. Skipping."
    cap log close
    exit 0
}

import delimited "`acs_csv'", clear ///
    stringcols(1) varnames(1)
foreach v in agg_earn n_earn ed_tot ed_ba pop {
    capture confirm numeric variable `v'
    if _rc destring `v', replace force
}
quietly count
local n_raw = r(N)
keep if !mi(agg_earn) & !mi(n_earn) & !mi(ed_tot) & !mi(ed_ba) & !mi(pop) ///
    & n_earn > 0 & ed_tot > 0 & agg_earn > 0
quietly count
local n_acs = r(N)
display as text "assess_09: ACS complete-case tracts: `n_acs' of `n_raw'."
tempfile acsT
save "`acsT'"

use geoid hc_cluster leiden_cluster ///
    using "$temp/final_assignments_${run_tag}_tract.dta", clear
merge 1:1 geoid using "`acsT'", keep(match) nogenerate
quietly count
local n_base = r(N)
display as text "assess_09: tracts in both clustering and ACS complete cases: `n_base'."

gen str5 county = substr(geoid, 1, 5)
gen str2 state  = substr(geoid, 1, 2)
gen byte nation = 1
gen str12 region = ""
foreach r in northeast midwest south west {
    foreach s of global states_`r' {
        quietly replace region = "`r'" if state == "`s'"
    }
}
quietly replace region = "west" if inlist(state, "02", "15") & region == ""
gen str12 division = ""
foreach d of global division_keys {
    foreach s of global div_`d' {
        quietly replace division = "`d'" if state == "`s'"
    }
}
tempfile base
save "`base'"

* County-level clusters, so that the ladder can be read at both scales from the
* same tract records.  These clusters are built from counties and therefore
* inherit county boundaries; see the discussion in Section A.14.
capture confirm file "$temp/final_assignments_${run_tag}_county.dta"
if _rc {
    display as text "assess_09: no county assignments; county-cluster rows skipped."
    use "`base'", clear
    gen long hc_county     = .
    gen long leiden_county = .
    save "`base'", replace
}
else {
    tempfile coT
    use geoid hc_cluster leiden_cluster ///
        using "$temp/final_assignments_${run_tag}_county.dta", clear
    rename geoid          county
    rename hc_cluster     hc_county
    rename leiden_cluster leiden_county
    save "`coT'"

    use "`base'", clear
    merge m:1 county using "`coT'", keep(master match) keepusing(hc_county leiden_county) nogenerate
    quietly count if mi(hc_county)
    if r(N) > 0 display as text "assess_09: `r(N)' tracts unmatched to county clusters."
    save "`base'", replace
}

tempname P
tempfile res
postfile `P' str16 key str16 levlabel str16 type long n_units ///
    double b double se double r2 double b_w double se_w ///
    using "`res'", replace

local keys "nation region division state hc_county leiden_county county hc leiden tract"
foreach k of local keys {
    if "`k'" == "nation"        local gv nation
    if "`k'" == "region"        local gv region
    if "`k'" == "division"      local gv division
    if "`k'" == "state"         local gv state
    if "`k'" == "hc_county"     local gv hc_county
    if "`k'" == "leiden_county" local gv leiden_county
    if "`k'" == "county"        local gv county
    if "`k'" == "hc"            local gv hc_cluster
    if "`k'" == "leiden"        local gv leiden_cluster
    if "`k'" == "tract"         local gv geoid

    if "`k'" == "nation"        local lab "National"
    if "`k'" == "region"        local lab "Census region"
    if "`k'" == "division"      local lab "Census division"
    if "`k'" == "state"         local lab "State"
    if "`k'" == "hc_county"     local lab "HC (county)"
    if "`k'" == "leiden_county" local lab "Leiden (county)"
    if "`k'" == "county"        local lab "County"
    if "`k'" == "hc"            local lab "HC (tract)"
    if "`k'" == "leiden"        local lab "Leiden (tract)"
    if "`k'" == "tract"         local lab "Tract"

    local ty "Administrative"
    if inlist("`k'", "hc", "leiden", "hc_county", "leiden_county") local ty "Functional"
    if "`k'" == "tract" local ty "Baseline"

    use "`base'", clear
    quietly drop if mi(`gv')
    collapse (sum) agg_earn n_earn ed_tot ed_ba pop, by(`gv')
    keep if n_earn > 0 & ed_tot > 0 & agg_earn > 0

    gen double y = agg_earn / n_earn
    gen double x = 100 * ed_ba / ed_tot

    quietly count
    if r(N) < 3 {
        display as text "  `lab' (`ty'): " r(N) " unit(s); slope undefined, posting missing."
        post `P' ("`k'") ("`lab'") ("`ty'") (r(N)) (.) (.) (.) (.) (.)
        continue
    }

    quietly regress y x
    local b  = _b[x]
    local se = _se[x]
    local n  = e(N)
    local r2 = e(r2)

    quietly regress y x [aweight = pop]
    local bw  = _b[x]
    local sew = _se[x]

    display as text "  `lab' (`ty'): units=" %7.0fc `n' ///
        "  b=" %7.1f `b' " (se " %5.1f `se' ")  R2=" %5.3f `r2'
    post `P' ("`k'") ("`lab'") ("`ty'") (`n') (`b') (`se') (`r2') (`bw') (`sew')
}
postclose `P'

use "`res'", clear
* rows are posted in ladder order (coarsest to finest); do not re-sort
gen lo   = b - invnormal(0.975) * se
gen hi   = b + invnormal(0.975) * se
gen lo_w = b_w - invnormal(0.975) * se_w
gen hi_w = b_w + invnormal(0.975) * se_w

label variable levlabel "Spatial unit"
label variable type     "Unit type"
label variable n_units  "Number of areal units (N)"
label variable b        "OLS slope: \$ mean earnings per +1pp BA"
label variable se       "Std. error (unweighted)"
label variable r2       "R-squared (unweighted)"
label variable b_w      "Slope, population-weighted"
format b se b_w se_w lo hi lo_w hi_w %9.1f
format r2 %5.3f

save           "$tables/maup_${run_tag}.dta", replace
export delimited using "$tables/maup_${run_tag}.csv", replace
display as result "  -> $tables/maup_${run_tag}.{dta,csv}"

display as text _n "{hline 72}"
display as text "MAUP: slope of mean earnings on % BA+ across spatial units"
display as text "{hline 72}"
list levlabel type n_units b se r2, noobs

if "${run_tag}" == "baseline" {
    capture file close _tex
    file open _tex using "$tables/maup.tex", write replace
    file write _tex "% Auto-generated by assess_09_maup.do; do not edit by hand." _n
    local dollar = char(92) + char(36)
    file write _tex "\begin{tabular}{llrrrr}" _n
    file write _tex "\toprule" _n
    file write _tex "             &      & Areal  & \multicolumn{2}{c}{Slope: `dollar' earnings / +1pp BA} & \\" _n
    file write _tex "Spatial unit & Type & units  & Unweighted & Pop.-weighted & \(R^2\) \\" _n
    file write _tex "\midrule" _n
    local _N = _N
    forvalues r = 1/`_N' {
        local lv  = levlabel[`r']
        local tp  = type[`r']
        local nu  : display %9.0fc n_units[`r']
        local nu = subinstr(trim("`nu'"), ",", "{,}", .)
        * a single areal unit admits no slope; print rules, not Stata's dot
        if mi(b[`r']) {
            file write _tex "`lv' & `tp' & `nu' & --- & --- & --- \\" _n
            continue
        }
        local bb  : display %9.0fc b[`r']
        local ss  : display %9.0fc se[`r']
        local bw  : display %9.0fc b_w[`r']
        local rr  : display %5.3f  r2[`r']
        local bb = subinstr(trim("`bb'"), ",", "{,}", .)
        local ss = subinstr(trim("`ss'"), ",", "{,}", .)
        local bw = subinstr(trim("`bw'"), ",", "{,}", .)
        local rr = trim("`rr'")
        file write _tex "`lv' & `tp' & `nu' & `bb' (`ss') & `bw' & `rr' \\" _n
    }
    file write _tex "\bottomrule" _n
    file write _tex "\end{tabular}" _n
    file close _tex
    display as result "  -> maup.tex (LaTeX fragment for tab:maup)"
}
else display as text "tex skipped: run_tag=${run_tag} != baseline"

display as result _n "assess_09 done."
cap log close
