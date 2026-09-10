clear all
set more off
cap log close
local _logdate = subinstr("$S_DATE", " ", "", .)
log using "$log/assess_12_random_null_${run_tag}_`_logdate'.log", replace text

* assess_12_random_null.do -- run the assess_11_robinson.do regressions on random partitions matched on the number of units, so the functional units can be judged against partitions that carry no commuting information.
* Called by: _master.do (assess phase, do_assess switch).
* Reads:  includes _inc/census_regions.do and _inc/robinson_engine.do; $temp/acs_pairs_tract_2023.csv; $temp/final_assignments_<run_tag>_tract.dta; $temp/random_nulls/random_<kind>_k<K>.dta (built by shell "$py" _py/random_partitions.py from $temp/tract_neighbors.dta when absent)
* Writes: $tables/robinson_null_<run_tag>.dta; $tables/robinson_null_<run_tag>.csv; via shell "$py" _py/robinson_figures.py --null, $figures/robinson_<name>.png including robinson_null.png
* Notes:  Two nulls are drawn: unconstrained random labels, which induce no aggregation bias, and contiguous regions grown on the tract adjacency graph, which keep spatial coherence but not commuting. The second is the one that isolates the role of boundaries. Defaults are $null_kinds, $null_ks (644 3100 4298 6627) and $null_reps (20); none is set in _master.do.
*
* assess_11_robinson.do shows that aggregating tracts into commuting-based
* units reproduces the tract-level estimate far better than aggregating them
* into counties at a matched number of units. The obvious reply is that any
* partition into several thousand pieces would do as well, and the commuting
* flows are doing no work. This script answers it by running the same
* regressions on partitions that carry no information, both matched on the
* number of units (random_partitions.py). The unconstrained null assigns tracts
* to clusters at random; because assignment is independent of everything, the
* tract-level estimate should survive aggregation almost intact. The contiguous
* null grows regions from random seeds on the tract adjacency graph, so it is
* spatially connected and compact like a real geography but blind to commuting.
* The comparison to draw is not "functional beats random". The question is
* where the county sits relative to the contiguous null: if the county is an
* outlier against partitions that share its spatial character but not its
* boundaries, then its boundaries, not its coarseness, distort the estimate.
* Estimation runs through _inc/robinson_engine.do, identical to assess_11, so
* the null and the estimates it benchmarks differ only in the partition.

include "$program/_inc/census_regions.do"

if "$null_kinds" == "" global null_kinds "unconstrained contiguous"
if "$null_ks"    == "" global null_ks    "644 3100 4298 6627"
if "$null_reps"  == "" global null_reps  20
local nulldir "$temp/random_nulls"

local acs_csv "$temp/acs_pairs_tract_2023.csv"
capture confirm file "`acs_csv'"
if _rc {
    display as error "assess_12: `acs_csv' missing; run assess_11 first. Skipping."
    cap log close
    exit 0
}
capture confirm file "`nulldir'/random_contiguous_k3100.dta"
if _rc {
    display as text "assess_12: generating random partitions ..."
    shell "$py" "$program/_py/random_partitions.py" ///
        --neighbors "$temp/tract_neighbors.dta" ///
        --universe  "$temp/final_assignments_${run_tag}_tract.dta" ///
        --outdir    "`nulldir'" --reps $null_reps --seed $seed
}
capture confirm file "`nulldir'/random_contiguous_k3100.dta"
if _rc {
    display as error "assess_12: random partitions unavailable; skipping."
    cap log close
    exit 0
}

* ------------------------------------------------------------ analysis data --
local sums "pop foreign ed_tot ed_lths ed_ba workers transit pov_univ pov_below hu vacant occ renter agg_rent clf unemp agg_hhinc hh agg_earn n_earn aland"

import delimited "`acs_csv'", clear stringcols(1) varnames(1)
foreach v of local sums {
    capture confirm numeric variable `v'
    if _rc destring `v', replace force
}
* complete cases are taken per relationship (_robsample), as in assess_11
gen double aland_km2 = aland / 1e6
drop aland
local sums = subinstr("`sums'", "aland", "aland_km2", .)
tempfile acsT
save "`acsT'"

use geoid using "$temp/final_assignments_${run_tag}_tract.dta", clear
merge 1:1 geoid using "`acsT'", keep(match) nogenerate
gen str2 state = substr(geoid, 1, 2)
encode state, gen(stnum)
quietly count
display as text "assess_12: analysis tracts: " r(N) "."
tempfile base
save "`base'"

include "$program/_inc/robinson_engine.do"

* ------------------------------------------------------------------- nulls --
tempname P
tempfile res
postfile `P' str16 rel str16 nullkind long k_target byte rep ///
    long n_units long nclust double b double se ///
    using "`res'", replace

foreach kind of global null_kinds {
    foreach K of global null_ks {
        local f "`nulldir'/random_`kind'_k`K'.dta"
        capture confirm file "`f'"
        if _rc {
            display as text "  `kind' k=`K': file absent; skipped."
            continue
        }
        forvalues r = 1/$null_reps {
            use "`base'", clear
            capture merge 1:1 geoid using "`f'", keepusing(c`r') keep(match) nogenerate
            if _rc {
                display as text "  `kind' k=`K' rep `r': merge failed; skipped."
                continue
            }
            quietly drop if mi(c`r')
            tempfile repT
            save "`repT'"

            foreach k of global rel_keys {
                use "`repT'", clear
                _robsample `k'
                collapse (sum) `sums' (firstnm) stnum, by(c`r')
                _robvars `k'
                _robfit y x
                post `P' ("`k'") ("`kind'") (`K') (`r') (r(n)) (r(nclust)) ///
                    (r(b)) (r(se))
            }
        }
        display as text "  `kind' k=`K': $null_reps replicates done."
    }
}
postclose `P'

use "`res'", clear
gen double t = b / se
gen str8 sign = "."
quietly replace sign = "null" if !mi(t) & abs(t) <  invnormal(0.975)
quietly replace sign = "pos"  if !mi(t) & abs(t) >= invnormal(0.975) & b > 0
quietly replace sign = "neg"  if !mi(t) & abs(t) >= invnormal(0.975) & b < 0
gen str64 rellab = ""
foreach k of global rel_keys {
    quietly replace rellab = "${l_`k'}" if rel == "`k'"
}
label variable k_target  "Target number of units"
label variable n_units   "Realised number of units"
label variable nullkind  "Random partition type"
save           "$tables/robinson_null_${run_tag}.dta", replace
export delimited using "$tables/robinson_null_${run_tag}.csv", replace
display as result "  -> $tables/robinson_null_${run_tag}.{dta,csv}"

capture confirm file "$tables/robinson_ladder_${run_tag}.csv"
if !_rc {
    display as text _n "assess_12: redrawing the Robinson figures with the null rows ..."
    shell "$py" "$program/_py/robinson_figures.py" ///
        --ladder "$tables/robinson_ladder_${run_tag}.csv" ///
        --specs  "$tables/robinson_specs_${run_tag}.csv" ///
        --null   "$tables/robinson_null_${run_tag}.csv" ///
        --outdir "$figures"
}

* ------------------------------------------------------------------ summary --
display as text _n "{hline 78}"
display as text "Null distribution of the slope, by relationship and partition type"
display as text "{hline 78}"
foreach k of global rel_keys {
    foreach kind of global null_kinds {
        quietly summarize b if rel == "`k'" & nullkind == "`kind'" & k_target == 3100
        if r(N) == 0 continue
        display as text "  " %-34s "${l_`k'}" " " %-14s "`kind'" ///
            "  mean=" %9.1f r(mean) "  sd=" %8.1f r(sd) ///
            "  [" %9.1f r(min) "," %9.1f r(max) "]"
    }
}

display as result _n "assess_12 done."
cap log close
