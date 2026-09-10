clear all
set more off
cap log close
local _logdate = subinstr("$S_DATE", " ", "", .)
log using "$log/assess_11_robinson_${run_tag}_`_logdate'.log", replace text

* assess_11_robinson.do -- test whether the areal unit changes the sign or significance of a screened set of ACS relationships, on the assess_09_maup.do ladder and then on every delineation produced.
* Called by: _master.do (assess phase, do_assess switch).
* Reads:  includes _inc/census_regions.do and _inc/robinson_engine.do; $temp/acs_pairs_tract_2023.csv (fetched by shell "$py" _py/fetch_acs_pairs.py when absent); $temp/final_assignments_<run_tag>_tract.dta; $temp/final_assignments_<run_tag>_county.dta; $temp/final_assignments_<tag>_tract.dta for every tag on disk
* Writes: $tables/robinson_ladder_<run_tag>.dta and .csv; $tables/robinson_specs_<run_tag>.dta and .csv; $tables/robinson_heatmap.tex (baseline run only); via shell "$py" _py/robinson_figures.py, $figures/robinson_<name>.png (assess_12 redraws them once the null rows exist)
* Notes:  Standard errors are clustered on state where a specification has enough states, because aggregated areal residuals are spatially correlated and the figures colour cells by significance. The specification pass ($rob_speccurve, default 1) sweeps every final_assignments_*_tract.dta in $temp, so any rerun left there enters the specification curve.
*
* assess_09_maup.do shows how one relationship, mean earnings on percent BA,
* shifts in magnitude across ten definitions of the areal unit; its estimate
* never changes sign. This script asks the sharper question: across a screened
* set of relationships, does the choice of unit change the sign or the
* significance of the estimate, and so the conclusion rather than the number?
* The design is Robinson's (1950). One set of tract records is aggregated many
* ways and the same regression is run on each, so every difference between
* rows is attributable to the unit alone. The ladder includes units this
* pipeline built from commuting flows, at granularities matched to the
* administrative units they are compared with. Holding the number of units
* fixed, any remaining difference is zonation, not scale.
* Every ACS quantity is an aggregate or a count, never a published median or
* rate, because only sums may be aggregated to arbitrary units; rates and means
* are formed after aggregation. We verified the table line numbers in
* fetch_acs_pairs.py against published national 2019-2023 totals before use:
* foreign born 13.76%, public transport to work 3.49%, poverty 12.73%, vacant
* units 10.55%, renter-occupied 34.95%, unemployment 5.25%, less than high
* school 10.70%, BA or higher 34.94%, and households 128,702,462, matching
* occupied units in B25003 exactly.
* Specifications with fewer than $rob_minn units are estimated but flagged: a
* significance claim resting on four Census regions is not one we ask a reader
* to accept.

include "$program/_inc/census_regions.do"

* Units with fewer than this many areas are reported but marked low-power.
if "$rob_minn" == "" global rob_minn 30

* ---------------------------------------------------------------- ACS data --
local acs_csv "$temp/acs_pairs_tract_2023.csv"
capture confirm file "`acs_csv'"
if _rc {
    display as text "assess_11: fetching ACS pair data ..."
    shell "$py" "$program/_py/fetch_acs_pairs.py" "`acs_csv'"
}
capture confirm file "`acs_csv'"
if _rc {
    display as error "assess_11: `acs_csv' unavailable; skipping."
    cap log close
    exit 0
}

capture confirm file "$temp/final_assignments_${run_tag}_tract.dta"
if _rc {
    display as error "assess_11: tract assignments missing; skipping."
    cap log close
    exit 0
}

* Summable ACS aggregates.  Nothing here is a rate.
local sums "pop foreign ed_tot ed_lths ed_ba workers transit pov_univ pov_below hu vacant occ renter agg_rent clf unemp agg_hhinc hh agg_earn n_earn aland"

import delimited "`acs_csv'", clear stringcols(1) varnames(1)
foreach v of local sums {
    capture confirm numeric variable `v'
    if _rc destring `v', replace force
}
quietly count
display as text "assess_11: ACS tracts read: " r(N) ///
    ". Complete cases are taken per relationship (_robsample)."
* land area in square kilometres, so density is people per km2
gen double aland_km2 = aland / 1e6
drop aland
local sums = subinstr("`sums'", "aland", "aland_km2", .)
tempfile acsT
save "`acsT'"

* ------------------------------------------------------- unit definitions --
use geoid hc_cluster leiden_cluster ///
    using "$temp/final_assignments_${run_tag}_tract.dta", clear
merge 1:1 geoid using "`acsT'", keep(match) nogenerate
quietly count
display as text "assess_11: tracts in both clustering and ACS: " r(N) "."

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
* numeric state key for clustering
encode state, gen(stnum)

* county-scale functional units, built by grouping whole counties
capture confirm file "$temp/final_assignments_${run_tag}_county.dta"
if _rc {
    display as text "assess_11: no county assignments; county-cluster rows skipped."
    gen long hc_county = .
    gen long leiden_county = .
}
else {
    tempfile coT
    preserve
    use geoid hc_cluster leiden_cluster ///
        using "$temp/final_assignments_${run_tag}_county.dta", clear
    rename geoid      county
    rename hc_cluster hc_county
    rename leiden_cluster leiden_county
    save "`coT'"
    restore
    merge m:1 county using "`coT'", keep(master match) ///
        keepusing(hc_county leiden_county) nogenerate
}
tempfile base
save "`base'"

include "$program/_inc/robinson_engine.do"

* ------------------------------------------------- pass 1: the unit ladder --
tempname P
tempfile res
postfile `P' str16 rel str24 unit str16 family long n_units long nclust ///
    double b double se double r2 double b_w double se_w ///
    using "`res'", replace

local units "nation region division state hc_county leiden_county county hc leiden tract"
foreach u of local units {
    if "`u'" == "nation"        local gv nation
    if "`u'" == "region"        local gv region
    if "`u'" == "division"      local gv division
    if "`u'" == "state"         local gv state
    if "`u'" == "hc_county"     local gv hc_county
    if "`u'" == "leiden_county" local gv leiden_county
    if "`u'" == "county"        local gv county
    if "`u'" == "hc"            local gv hc_cluster
    if "`u'" == "leiden"        local gv leiden_cluster
    if "`u'" == "tract"         local gv geoid

    if "`u'" == "nation"        local ulab "National"
    if "`u'" == "region"        local ulab "Census region"
    if "`u'" == "division"      local ulab "Census division"
    if "`u'" == "state"         local ulab "State"
    if "`u'" == "hc_county"     local ulab "HC (county)"
    if "`u'" == "leiden_county" local ulab "Leiden (county)"
    if "`u'" == "county"        local ulab "County"
    if "`u'" == "hc"            local ulab "HC (tract)"
    if "`u'" == "leiden"        local ulab "Leiden (tract)"
    if "`u'" == "tract"         local ulab "Tract"

    local fam "Administrative"
    if inlist("`u'", "hc", "leiden", "hc_county", "leiden_county") local fam "Functional"
    if "`u'" == "tract" local fam "Baseline"

    use "`base'", clear
    quietly drop if mi(`gv')
    quietly count
    if r(N) == 0 {
        display as text "  `ulab': no rows; skipped."
        continue
    }
    tempfile unitT
    save "`unitT'"

    * one collapse per relationship: each is estimated on its own complete
    * cases (_robsample), so the sample matches assess_09/assess_10
    * relationship for relationship.
    foreach k of global rel_keys {
        use "`unitT'", clear
        _robsample `k'
        collapse (sum) `sums' (firstnm) stnum, by(`gv')
        _robvars `k'
        _robfit y x
        local n = r(n)
        local b = r(b)
        local se = r(se)
        local r2 = r(r2)
        local nc = r(nclust)
        _robfit y x, weightvar(pop)
        local bw = r(b)
        local sew = r(se)
        post `P' ("`k'") ("`ulab'") ("`fam'") (`n') (`nc') ///
            (`b') (`se') (`r2') (`bw') (`sew')
        if "`k'" == "earn_ba" local n_first = `n'
    }
    display as text "  `ulab': done (`n_first' units on earn_ba)."
}
postclose `P'

use "`res'", clear
gen double t = b / se
gen str8 sign = "."
quietly replace sign = "null" if !mi(t) & abs(t) <  invnormal(0.975)
quietly replace sign = "pos"  if !mi(t) & abs(t) >= invnormal(0.975) & b > 0
quietly replace sign = "neg"  if !mi(t) & abs(t) >= invnormal(0.975) & b < 0
gen byte lowpower = n_units < $rob_minn
gen str64 rellab = ""
foreach k of global rel_keys {
    quietly replace rellab = "${l_`k'}" if rel == "`k'"
}
label variable n_units  "Number of areal units"
label variable nclust   "States available for clustering"
label variable b        "Slope (state-clustered s.e.)"
label variable sign     "pos / null / neg at 5%"
label variable lowpower "Fewer than $rob_minn units"
save           "$tables/robinson_ladder_${run_tag}.dta", replace
export delimited using "$tables/robinson_ladder_${run_tag}.csv", replace
display as result "  -> $tables/robinson_ladder_${run_tag}.{dta,csv}"

display as text _n "{hline 78}"
display as text "Sign of the estimate by relationship and unit (5%, state-clustered)"
display as text "{hline 78}"
list rellab unit n_units b se sign if rel == "earn_rent", noobs

* ------------------------------------ pass 2: every delineation the pipeline
* produced.  The unit ladder above compares one delineation against the
* administrative maps; this pass re-runs it on all of them, so that the
* specification curve reports the full multiverse rather than a chosen point in
* it.  Robustness checks, calibration perturbations, temporal windows and
* threshold variants all enter on the same footing.
if "$rob_speccurve" == "" global rob_speccurve 1

if $rob_speccurve {
    local tagfiles : dir "$temp" files "final_assignments_*_tract.dta"
    local ntag : word count `tagfiles'
    display as text _n "assess_11: specification pass over `ntag' delineations."

    tempname Q
    tempfile spec
    postfile `Q' str16 rel str32 tag str8 method long n_units long nclust ///
        double b double se double b_w double se_w ///
        using "`spec'", replace

    foreach f of local tagfiles {
        local tag = subinstr("`f'", "final_assignments_", "", .)
        local tag = subinstr("`tag'", "_tract.dta", "", .)

        capture noisily {
            preserve
            use geoid hc_cluster leiden_cluster using "$temp/`f'", clear
            tempfile aT
            save "`aT'"
            restore
        }
        if _rc {
            display as text "  `tag': unreadable; skipped."
            continue
        }

        foreach m in hc leiden {
            if "`m'" == "hc"     local gv hc_cluster
            if "`m'" == "leiden" local gv leiden_cluster

            use "`base'", clear
            capture drop hc_cluster leiden_cluster
            merge 1:1 geoid using "`aT'", keep(match) nogenerate
            quietly drop if mi(`gv')
            quietly count
            if r(N) == 0 continue
            tempfile specT
            save "`specT'"

            foreach k of global rel_keys {
                use "`specT'", clear
                _robsample `k'
                collapse (sum) `sums' (firstnm) stnum, by(`gv')
                _robvars `k'
                _robfit y x
                local n = r(n)
                local b = r(b)
                local se = r(se)
                local nc = r(nclust)
                _robfit y x, weightvar(pop)
                post `Q' ("`k'") ("`tag'") ("`m'") (`n') (`nc') ///
                    (`b') (`se') (r(b)) (r(se))
            }
        }
        display as text "  `tag': done."
    }
    postclose `Q'

    use "`spec'", clear
    gen double t = b / se
    gen str8 sign = "."
    quietly replace sign = "null" if !mi(t) & abs(t) <  invnormal(0.975)
    quietly replace sign = "pos"  if !mi(t) & abs(t) >= invnormal(0.975) & b > 0
    quietly replace sign = "neg"  if !mi(t) & abs(t) >= invnormal(0.975) & b < 0
    gen str64 rellab = ""
    foreach k of global rel_keys {
        quietly replace rellab = "${l_`k'}" if rel == "`k'"
    }
    save           "$tables/robinson_specs_${run_tag}.dta", replace
    export delimited using "$tables/robinson_specs_${run_tag}.csv", replace
    display as result "  -> $tables/robinson_specs_${run_tag}.{dta,csv}"

    display as text _n "{hline 78}"
    display as text "Delineations per relationship, and how many change sign"
    display as text "{hline 78}"
    quietly levelsof rel, local(_rels)
    foreach k of local _rels {
        quietly count if rel == "`k'" & sign == "pos"
        local np = r(N)
        quietly count if rel == "`k'" & sign == "neg"
        local nn = r(N)
        quietly count if rel == "`k'" & sign == "null"
        local n0 = r(N)
        display as text "  " %-38s "${l_`k'}" ///
            "  pos=" %4.0f `np' "  null=" %4.0f `n0' "  neg=" %4.0f `nn'
    }
}
else display as text "assess_11: specification pass disabled (rob_speccurve=0)."

* ------------------------------------------------------------ LaTeX output --
if "${run_tag}" == "baseline" {
    use "$tables/robinson_ladder_${run_tag}.dta", clear
    * paper table: one row per relationship, one column per unit, sign only
    local u1 "Census region"
    local u2 "Census division"
    local u3 "State"
    local u4 "HC (county)"
    local u5 "Leiden (county)"
    local u6 "County"
    local u7 "HC (tract)"
    local u8 "Leiden (tract)"
    local u9 "Tract"
    * A LaTeX row ends with two backslashes.  Build every backslash from
    * char(92) so this block does not depend on how backslashes survive
    * Stata string literals and macro expansion.
    local bs = char(92)
    local eol "`bs'`bs'"
    capture file close _tex
    file open _tex using "$tables/robinson_heatmap.tex", write replace
    file write _tex "% Auto-generated by assess_11_robinson.do; do not edit by hand." _n
    file write _tex "% `bs'posmark, `bs'nullmark and `bs'negmark are defined in appendixstyle.sty." _n
    file write _tex "`bs'begin{tabular}{l*{9}{c}}" _n
    file write _tex "`bs'toprule" _n
    file write _tex " & `bs'multicolumn{3}{c}{Administrative, coarse} & `bs'multicolumn{2}{c}{Functional (county)} & & `bs'multicolumn{2}{c}{Functional (tract)} & `eol'" _n
    file write _tex "`bs'cmidrule(lr){2-4}`bs'cmidrule(lr){5-6}`bs'cmidrule(lr){8-9}" _n
    file write _tex "Relationship & Reg & Div & St & HC & Leiden & County & HC & Leiden & Tract `eol'" _n
    file write _tex "`bs'midrule" _n
    foreach k of global rel_keys {
        * "~" in the stored label is a tilde, not a non-breaking space
        local line = subinstr("${l_`k'}", "~", "`bs'ensuremath{`bs'sim}", .)
        forvalues c = 1/9 {
            local uu "`u`c''"
            quietly levelsof sign if rel == "`k'" & unit == "`uu'", local(_s) clean
            local cell "`bs'ensuremath{`bs'cdot}"
            if "`_s'" == "pos"  local cell "`bs'posmark"
            if "`_s'" == "null" local cell "`bs'nullmark"
            if "`_s'" == "neg"  local cell "`bs'negmark"
            local line "`line' & `cell'"
        }
        file write _tex "`line' `eol'" _n
    }
    file write _tex "`bs'bottomrule" _n
    file write _tex "`bs'end{tabular}" _n
    file close _tex
    display as result "  -> robinson_heatmap.tex"
}

* ---------------------------------------------------------------- figures --
capture confirm file "$tables/robinson_ladder_${run_tag}.csv"
if !_rc {
    display as text _n "assess_11: drawing figures ..."
    shell "$py" "$program/_py/robinson_figures.py" ///
        --ladder "$tables/robinson_ladder_${run_tag}.csv" ///
        --specs  "$tables/robinson_specs_${run_tag}.csv" ///
        --outdir "$figures"
}

display as result _n "assess_11 done."
cap log close
