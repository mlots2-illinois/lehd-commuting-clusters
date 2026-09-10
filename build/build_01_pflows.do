clear all
cap log close
local _logdate = subinstr("$S_DATE", " ", "", .)
log using "$log/build_01_pflows_${run_tag}_`_logdate'.log", replace text

* build_01_pflows.do -- converts each pooled OD file into an undirected pair file with the Tolbert-Killian proportional flow P_ij = f_ij / min(rLF_i, rLF_j).
* Called by: _master.do (build phase; once per job type in $job_types with $pool_tag set, then once per $min_flow cut in 2 3 5 10 20 with $flow_tag = _mf<m>).
* Reads:  $temp/5yr/od_s000_<lvl>_5yr_<ys>_<ye><pool_tag>.dta; $temp/3yr/od_s000_<lvl>_3yr_<ys>_<ye><pool_tag>.dta; $temp/1yr/od_s000_<lvl>_1yr_<year><pool_tag>.dta
* Writes: $temp/pflow_<lvl>_5yr_<ys>_<ye><pool_tag><flow_tag>.dta; $temp/pflow_<lvl>_3yr_<ys>_<ye><pool_tag><flow_tag>.dta; $temp/pflow_<lvl>_1yr_<year><pool_tag><flow_tag>.dta
* Notes:  rLF, the resident labor force, is the sum of S000 over all work destinations of a home unit and is computed before within-unit flows are dropped. AK, HI and the territories in $noncontig_fips leave at this step. We drop pairs with f_ij < $min_flow, so the _mf<m> files are thresholded copies of the baseline; with $min_flow = 0 only zero-flow pairs go.

capture program drop compute_rLF
program define compute_rLF
    syntax , geo(string)
    if !inlist("`geo'", "county", "tract") {
        di as error "geo() must be 'county' or 'tract'"
        exit 198
    }
    local hvar "h_`geo'"
    local wvar "w_`geo'"
    confirm variable `hvar' `wvar' s000

    foreach _fp of global noncontig_fips {
        drop if substr(`hvar',1,2) == "`_fp'" | substr(`wvar',1,2) == "`_fp'"
    }

    capture drop rLF
    gegen double rLF = total(s000), by(`hvar')
    label variable rLF "Resident labor force (sum S000 by `hvar')"
end

capture program drop compute_pflows
program define compute_pflows
    syntax , geo(string) saveas(string)
    if !inlist("`geo'", "county", "tract") {
        di as error "geo() must be 'county' or 'tract'"
        exit 198
    }
    local hvar "h_`geo'"
    local wvar "w_`geo'"
    confirm variable `hvar' `wvar' s000 rLF

    preserve
    bysort `hvar': keep if _n == 1
    keep `hvar' rLF
    rename `hvar' geo_fips
    rename rLF rLF_geo
    tempfile rlf_lookup
    save `rlf_lookup'
    restore

    drop if `hvar' == `wvar'
    gen str geo_i = cond(`hvar' < `wvar', `hvar', `wvar')
    gen str geo_j = cond(`hvar' < `wvar', `wvar', `hvar')

    gcollapse (sum) f_ij = s000, by(geo_i geo_j)

    rename geo_i geo_fips
    merge m:1 geo_fips using `rlf_lookup', keep(master match)
    quietly count if _merge == 1
    local _n_unmatched = r(N)
    display as text "  pairs dropped for missing rLF (geo_i endpoint): `_n_unmatched' of " _N
    if `_n_unmatched' > 0.001 * _N {
        display as error "WARNING: >0.1% of pairs lack an rLF match on geo_i; check universe alignment (pool vs rLF)."
    }
    drop if _merge == 1
    drop _merge
    rename rLF_geo rLF_i
    rename geo_fips geo_i

    rename geo_j geo_fips
    merge m:1 geo_fips using `rlf_lookup', keep(master match)
    quietly count if _merge == 1
    local _n_unmatched = r(N)
    display as text "  pairs dropped for missing rLF (geo_j endpoint): `_n_unmatched' of " _N
    if `_n_unmatched' > 0.001 * _N {
        display as error "WARNING: >0.1% of pairs lack an rLF match on geo_j; check universe alignment (pool vs rLF)."
    }
    drop if _merge == 1
    drop _merge
    rename rLF_geo rLF_j
    rename geo_fips geo_j

    gen double min_rLF = min(rLF_i, rLF_j)
    gen double P_ij = f_ij / min_rLF if min_rLF > 0 & !mi(min_rLF)
    label variable P_ij "Proportional commuting flow (Tolbert-Killian)"
    label variable f_ij "Bidirectional flow (i->j + j->i)"

    drop if mi(P_ij) | P_ij == 0 | f_ij < ${min_flow}
    drop min_rLF
    order geo_i geo_j f_ij rLF_i rLF_j P_ij
    sort geo_i geo_j

    compress
    save "`saveas'", replace
    di as text "  Saved: `saveas' (" _N " pairs)"
end

local y3a = $yr_end - 2

foreach geo of global levels {

    di as text _n "=== `geo' / 5yr_${yr_start}_${yr_end}${pool_tag}${flow_tag} ==="
    use "$temp/5yr/od_s000_`geo'_5yr_${yr_start}_${yr_end}${pool_tag}.dta", clear
    compute_rLF, geo(`geo')
    compute_pflows, geo(`geo') ///
        saveas("$temp/pflow_`geo'_5yr_${yr_start}_${yr_end}${pool_tag}${flow_tag}.dta")

    di as text _n "=== `geo' / 3yr_`y3a'_${yr_end}${pool_tag}${flow_tag} ==="
    use "$temp/3yr/od_s000_`geo'_3yr_`y3a'_${yr_end}${pool_tag}.dta", clear
    compute_rLF, geo(`geo')
    compute_pflows, geo(`geo') ///
        saveas("$temp/pflow_`geo'_3yr_`y3a'_${yr_end}${pool_tag}${flow_tag}.dta")

    forvalues y = $yr_start/$yr_end {
        di as text _n "=== `geo' / 1yr_`y'${pool_tag}${flow_tag} ==="
        use "$temp/1yr/od_s000_`geo'_1yr_`y'${pool_tag}.dta", clear
        compute_rLF, geo(`geo')
        compute_pflows, geo(`geo') ///
            saveas("$temp/pflow_`geo'_1yr_`y'${pool_tag}${flow_tag}.dta")
    }
}

display as result _n "build_01 done."

cap log close
