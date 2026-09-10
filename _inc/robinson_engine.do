* robinson_engine.do -- shared estimation engine for the unit-of-analysis (MAUP) work: the relationship specifications, the sample rule and the regression that assess_11 and assess_12 both call.
* Included by: assess/assess_11_robinson.do (real delineations), assess/assess_12_random_null.do (random partitions).
* Expects: the ACS tract aggregates in memory (agg_earn n_earn ed_ba ed_tot pov_below pov_univ ed_lths foreign pop aland_km2 vacant hu renter occ agg_hhinc hh unemp clf) and stnum for clustering.
* Notes:  the null and the estimates it benchmarks run through the same code path; with separate implementations a difference between them could be an artefact of the estimator rather than of the partition. Defines $rel_keys, $r_<key>, $l_<key> and the programs _robfit, _robsample, _robvars.

* ------------------------------------------------------------ relationships --
* key -> ynum ydenom yscale xnum xdenom xscale.  A scale of LOG means the
* construct is the natural log of the ratio; otherwise it multiplies it.
global rel_keys "earn_ba earn_pov lths_fb dens_ba pov_vac earn_rent hhinc_rent earn_unemp hhinc_unemp ba_unemp"

global r_earn_ba     "agg_earn n_earn 1 ed_ba ed_tot 100"
global r_earn_pov    "agg_earn n_earn 1 pov_below pov_univ 100"
global r_lths_fb     "ed_lths ed_tot 100 foreign pop 100"
* Density is logged because tract densities span several orders of magnitude; the other constructs are shares (x100) or per-capita means (x1).
global r_dens_ba     "pop aland_km2 LOG ed_ba ed_tot 100"
global r_pov_vac     "pov_below pov_univ 100 vacant hu 100"
global r_earn_rent   "agg_earn n_earn 1 renter occ 100"
global r_hhinc_rent  "agg_hhinc hh 1 renter occ 100"
global r_earn_unemp  "agg_earn n_earn 1 unemp clf 100"
global r_hhinc_unemp "agg_hhinc hh 1 unemp clf 100"
global r_ba_unemp    "ed_ba ed_tot 100 unemp clf 100"

global l_earn_ba     "Mean earnings ~ pct BA+"
global l_earn_pov    "Mean earnings ~ poverty rate"
global l_lths_fb     "Pct < HS ~ pct foreign-born"
global l_dens_ba     "Log density ~ pct BA+"
global l_pov_vac     "Poverty rate ~ vacancy rate"
global l_earn_rent   "Mean earnings ~ pct renter-occupied"
global l_hhinc_rent  "Mean HH income ~ pct renter-occupied"
global l_earn_unemp  "Mean earnings ~ unemployment rate"
global l_hhinc_unemp "Mean HH income ~ unemployment rate"
global l_ba_unemp    "Pct BA+ ~ unemployment rate"

* ------------------------------------------------------------------ engine --
* Estimates one relationship on the data already collapsed to a unit, and
* posts the result.  Expects y and x to exist; clusters on stnum when the
* specification has enough states to support it.
capture program drop _robfit
program define _robfit, rclass
    syntax varlist(min=2 max=2) [, Weightvar(varname) ]
    tokenize `varlist'
    local y `1'
    local x `2'
    quietly count if !mi(`y') & !mi(`x')
    return scalar n = r(N)
    if r(N) < 3 {
        return scalar b = .
        return scalar se = .
        return scalar r2 = .
        return scalar nclust = 0
        exit
    }
    quietly levelsof stnum if !mi(`y') & !mi(`x'), local(_st)
    local nc : word count `_st'
    return scalar nclust = `nc'
    * Fewer than ten state clusters (regions, divisions) leave cluster-robust SEs with too few degrees of freedom, so we fall back to heteroskedasticity-robust SEs.
    if `nc' >= 10 {
        if "`weightvar'" == "" quietly regress `y' `x', vce(cluster stnum)
        else quietly regress `y' `x' [aweight = `weightvar'], vce(cluster stnum)
    }
    else {
        if "`weightvar'" == "" quietly regress `y' `x', robust
        else quietly regress `y' `x' [aweight = `weightvar'], robust
    }
    return scalar b  = _b[`x']
    return scalar se = _se[`x']
    return scalar r2 = e(r2)
end

* Restricts the TRACT file to relationship `k''s own complete cases: all four
* inputs reported and both denominators positive (numerators too where the
* construct is a log).  Applied before collapsing, so each relationship is
* estimated on every tract that can contribute to it -- the same rule that
* assess_09/assess_10 apply to the single earnings-on-education relationship,
* so the ladder, the parameter sweep and the null share one sample per
* relationship.
capture program drop _robsample
program define _robsample
    args k
    local p "${r_`k'}"
    local yn : word 1 of `p'
    local yd : word 2 of `p'
    local ys : word 3 of `p'
    local xn : word 4 of `p'
    local xd : word 5 of `p'
    local xs : word 6 of `p'
    quietly keep if !mi(`yn') & !mi(`yd') & !mi(`xn') & !mi(`xd') ///
        & `yd' > 0 & `xd' > 0
    * A zero numerator has no log; dropping it here rather than in _robvars keeps the sample identical across units for the same relationship.
    if "`ys'" == "LOG" quietly keep if `yn' > 0
    if "`xs'" == "LOG" quietly keep if `xn' > 0
end

* Builds y and x on the collapsed data for relationship `k'.
capture program drop _robvars
program define _robvars
    args k
    local p "${r_`k'}"
    local yn : word 1 of `p'
    local yd : word 2 of `p'
    local ys : word 3 of `p'
    local xn : word 4 of `p'
    local xd : word 5 of `p'
    local xs : word 6 of `p'
    capture drop y
    capture drop x
    if "`ys'" == "LOG" gen double y = log(`yn' / `yd') if `yd' > 0 & `yn' > 0
    else               gen double y = `ys' * `yn' / `yd' if `yd' > 0
    if "`xs'" == "LOG" gen double x = log(`xn' / `xd') if `xd' > 0 & `xn' > 0
    else               gen double x = `xs' * `xn' / `xd' if `xd' > 0
end
