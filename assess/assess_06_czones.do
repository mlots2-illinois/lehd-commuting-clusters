clear all
cap log close
local _logdate = subinstr("$S_DATE", " ", "", .)
log using "$log/assess_06_czones_${run_tag}_`_logdate'.log", replace text

* assess_06_czones.do -- compare each partition with the USDA/ERS commuting zones (CZ), the standard county-based delineation, by NMI, ARI and per-CZ purity.
* Called by: _master.do (assess phase, do_assess switch).
* Reads:  includes _inc/compare_partitions.do and _inc/maxshare_purity.do; $proj/03_data/01_raw/czones/county_to_cz.dta (or .csv); $temp/final_assignments_<run_tag>_<lvl>.dta
* Writes: $tables/cz_comparison_<run_tag>.dta; $tables/cz_comparison_<run_tag>.csv; $tables/cz_comparison.tex (baseline run only)
* Notes:  Commuting zones are defined on counties, so tracts inherit the zone of their county through the first five digits of the geoid. The crosswalk is not part of the bundle; when it is absent the script exits cleanly and the CZ table is simply not produced.

include "$program/_inc/compare_partitions.do"
include "$program/_inc/maxshare_purity.do"

local cz_dir "$proj/03_data/01_raw/czones"
local cz_dta "`cz_dir'/county_to_cz.dta"
local cz_csv "`cz_dir'/county_to_cz.csv"

tempfile czxw
capture confirm file "`cz_dta'"
if !_rc {
    use "`cz_dta'", clear
}
else {
    capture confirm file "`cz_csv'"
    if !_rc {
        import delimited using "`cz_csv'", varnames(1) stringcols(_all) clear
    }
    else {
        display as error "{hline 70}"
        display as error "assess_06: county->CZ crosswalk not found. Skipping CZ comparison."
        display as error "  Place a USDA/ERS commuting-zone crosswalk at:"
        display as error "    `cz_dta'   (or .csv)"
        display as error "  columns: geoid (5-digit county FIPS), cz (zone id)."
        display as error "{hline 70}"
        cap log close
        exit 0
    }
}
keep geoid cz
duplicates drop geoid, force
egen long cz_id = group(cz)
assert !missing(cz_id)
keep geoid cz_id
save `czxw'
quietly count
display as text "County->CZ crosswalk: " r(N) " counties."

tempfile _acc
local _have 0

foreach lvl of global levels {

    display as text _n "{hline 70}"
    display as text "===== assess_06 LEVEL: `lvl' ====="
    display as text "{hline 70}"

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

        if "`lvl'" == "tract" gen str5 _cty = substr(geoid, 1, 5)
        else                  gen str5 _cty = geoid
        rename geoid _unit
        rename _cty geoid
        quietly count
        local n_pre = r(N)
        merge m:1 geoid using `czxw', keep(match) keepusing(cz_id) nogenerate

        quietly count
        local n_units = r(N)
        display as text "  CZ merge (`lvl' `mlab'): `n_units' of `n_pre' units matched (" ///
            %5.1f 100 * `n_units' / `n_pre' "%)"
        if `n_units' <= 3000 {
            display as error "assess_06: only `n_units' units matched the CZ crosswalk;" ///
                " check the crosswalk vintage (county FIPS in county_to_cz must match" ///
                " the LODES county definitions, e.g. 2010 vs 2020 delineations)."
        }
        assert `n_units' > 3000
        quietly levelsof cz_id, local(_czs)
        local n_cz : word count `_czs'

        compare_partitions, a(clus) b(cz_id)
        local NMI = r(NMI)
        local ARI = r(ARI)

        maxshare_purity, refvar(cz_id) clustervar(clus) gen(purity)
        gen byte _full = (purity == 1)
        quietly summarize purity
        local mean_pur = r(mean)
        quietly summarize _full
        local pct_full = r(mean) * 100

        display as text "  `mlab' (`lvl'): n_cz=`n_cz'  NMI=" %5.3f `NMI' ///
            "  ARI=" %5.3f `ARI' "  mean purity=" %5.3f `mean_pur' ///
            "  %fully=" %4.1f `pct_full'

        clear
        set obs 1
        gen str6   geo_level   = "`lvl'"
        gen str8   method      = "`mlab'"
        gen long   n_cz        = `n_cz'
        gen long   n_units     = `n_units'
        gen double NMI         = `NMI'
        gen double ARI         = `ARI'
        gen double mean_purity = `mean_pur'
        gen double pct_fully   = `pct_full'

        if `_have' append using "`_acc'"
        save "`_acc'", replace
        local _have 1
    }
}

capture confirm file "`_acc'"
if _rc {
    display as error "No CZ comparison produced (missing assignments)."
    cap log close
    exit 0
}

use "`_acc'", clear
gen byte _ordg = cond(geo_level == "tract", 1, 2)
gen byte _ordm = cond(method == "HC", 1, 2)
sort _ordg _ordm
drop _ordg _ordm

label variable geo_level   "Geography (tract or county)"
label variable method      "Algorithm"
label variable n_cz        "Commuting zones covered"
label variable NMI         "Normalized MI with CZ partition"
label variable ARI         "Adjusted Rand index with CZ partition"
label variable mean_purity "Mean per-CZ max-share purity"
label variable pct_fully   "% CZs wholly within one cluster"
format NMI ARI mean_purity %5.3f
format pct_fully %5.1f

save           "$tables/cz_comparison_${run_tag}.dta", replace
export delimited using "$tables/cz_comparison_${run_tag}.csv", replace
display as result "  -> $tables/cz_comparison_${run_tag}.{csv,dta}"

list geo_level method n_cz NMI ARI mean_purity pct_fully, noobs

if "${run_tag}" == "baseline" {
    capture file close _tex
    file open _tex using "$tables/cz_comparison.tex", write replace
    file write _tex "% Auto-generated by assess_06_czones.do; do not edit by hand." _n
    file write _tex "\begin{tabular}{llrrrr}" _n
    file write _tex "\toprule" _n
    file write _tex "Geography & Method & \(n_{\text{CZ}}\) & NMI & ARI & \% fully preserved \\" _n
    file write _tex "\midrule" _n
    local _N = _N
    forvalues r = 1/`_N' {
        local geo = strproper(geo_level[`r'])
        local mth = method[`r']
        local nz  = n_cz[`r']
        local nm  : display %5.3f NMI[`r']
        local ar  : display %5.3f ARI[`r']
        local pf  : display %4.1f pct_fully[`r']
        file write _tex "`geo' & `mth' & `nz' & `=trim("`nm'")' & `=trim("`ar'")' & `=trim("`pf'")' \\" _n
    }
    file write _tex "\bottomrule" _n
    file write _tex "\end{tabular}" _n
    file close _tex
    display as result "  -> cz_comparison.tex (LaTeX fragment for tab:cz_comparison)"
}
else display as text "tex skipped: run_tag=${run_tag} != baseline"

display as result _n "assess_06 done."

cap log close
