clear all
cap log close
local _logdate = subinstr("$S_DATE", " ", "", .)
log using "$log/assess_03_temporal_${run_tag}_`_logdate'.log", replace text

* assess_03_temporal.do -- compare the baseline partition with the 3-year and 1-year window reruns, and the five single-year partitions with each other, to measure temporal stability.
* Called by: _master.do (assess phase, do_assess switch).
* Reads:  includes _inc/compare_partitions.do; $temp/final_assignments_<tag>_<lvl>.dta for tag = baseline, temporal_3yr_2021_2023, temporal_1yr_<year>
* Writes: $tables/temporal_agreement.dta; $tables/temporal_yoy.dta (gen/tab_temporal_stability.do and gen/tab_temporal_yoy.do turn these into the manuscript tables)
* Notes:  Each pair is compared on the intersection of the two unit sets. We print how many units each side loses so that a low score cannot be an artefact of coverage. The window reruns come from robust_01_window_3yr.do and robust_02_window_1yr.do; a level with any rerun missing is skipped rather than reported partially.

include "$program/_inc/compare_partitions.do"

local periods    "baseline temporal_3yr_2021_2023 temporal_1yr_2023"
local labels     "5yr 3yr 1yr"
local methods    "hc leiden"

tempfile out
local first 1

foreach lvl of global levels {

    display as text _n "{hline 70}"
    display as text "===== assess_03 LEVEL: `lvl' ====="
    display as text "{hline 70}"

    local missing ""
    foreach p of local periods {
        capture confirm file "$temp/final_assignments_`p'_`lvl'.dta"
        if _rc local missing "`missing' `p'"
    }
    if "`missing'" != "" {
        display as error "Missing temporal runs (`lvl'): `missing'; skipping this level."
        continue
    }

    local n_periods : word count `periods'

    forvalues i = 1/`n_periods' {
        forvalues j = `=`i' + 1'/`n_periods' {

            local p1   : word `i' of `periods'
            local p2   : word `j' of `periods'
            local lbl1 : word `i' of `labels'
            local lbl2 : word `j' of `labels'

            use geoid hc_cluster leiden_cluster ///
                using "$temp/final_assignments_`p1'_`lvl'.dta", clear
            rename hc_cluster     hc_a
            rename leiden_cluster le_a
            quietly count
            local n_p1 = r(N)
            tempfile a
            save `a'

            use geoid hc_cluster leiden_cluster ///
                using "$temp/final_assignments_`p2'_`lvl'.dta", clear
            rename hc_cluster     hc_b
            rename leiden_cluster le_b
            quietly count
            local n_p2 = r(N)
            merge 1:1 geoid using `a', keep(match) nogenerate
            quietly count
            local n_both = r(N)
            local n_dropped   = `n_p2' - `n_both'
            local n_dropped_a = `n_p1' - `n_both'
            if `n_dropped' > 0 {
                display as text "  [`lvl' `lbl1'-`lbl2'] " `n_dropped' ///
                    " of " `n_p2' " `lbl2'-units absent from `lbl1' " ///
                    "(compared on the " `n_both' "-unit intersection)"
            }
            if `n_dropped_a' > 0 {
                display as text "  [`lvl' `lbl1'-`lbl2'] " `n_dropped_a' ///
                    " of " `n_p1' " `lbl1'-units absent from `lbl2' " ///
                    "(compared on the " `n_both' "-unit intersection)"
            }

            foreach m of local methods {
                local va = cond("`m'" == "hc", "hc_a", "le_a")
                local vb = cond("`m'" == "hc", "hc_b", "le_b")

                compare_partitions, a(`va') b(`vb')

                preserve
                clear
                quietly set obs 1
                gen str6   geo_level = "`lvl'"
                gen str8   method    = "`m'"
                gen str20  period_a  = "`lbl1'"
                gen str20  period_b  = "`lbl2'"
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
        }
    }
}

capture confirm file "`out'"
if _rc {
    display as error "No temporal comparisons produced; both levels missing prerequisites."
    display as error "  robust_01_window_3yr.do and robust_02_window_1yr.do produce them (do_robust_window in _master.do)."
    cap log close
    exit 0
}

use `out', clear
sort geo_level method period_a period_b

label variable geo_level "Geography (tract or county)"
label variable method    "Algorithm"
label variable period_a  "Earlier period"
label variable period_b  "Later period"
label variable n         "Units compared (intersection of both pools)"
label variable n_dropped "Units in later pool but absent from earlier pool"
label variable NMI       "Normalized MI (arithmetic; not chance-corrected)"
label variable ARI       "Adjusted Rand index (chance-corrected)"
label variable AMI       "Adjusted mutual information (chance-corrected)"

format NMI ARI AMI %7.4f

save "$tables/temporal_agreement.dta", replace

display as text _n "{hline 70}"
display as text "Temporal agreement (pairwise across pool windows)"
display as text "{hline 70}"
list geo_level method period_a period_b n NMI ARI, noobs sepby(geo_level method)

local yrs     "2019 2020 2021 2022 2023"
local methods "hc leiden"

tempfile yout
local yfirst 1

foreach lvl of global levels {

    display as text _n "{hline 70}"
    display as text "===== assess_03 YoY LEVEL: `lvl' ====="
    display as text "{hline 70}"

    local ymiss ""
    foreach y of local yrs {
        capture confirm file "$temp/final_assignments_temporal_1yr_`y'_`lvl'.dta"
        if _rc local ymiss "`ymiss' `y'"
    }
    if "`ymiss'" != "" {
        display as error "Missing single-year runs (`lvl'):`ymiss'; skipping YoY for this level."
        display as error "  Produce them via robust_02_window_1yr.do (set \$do_robust = 1 and"
        display as error "  \$do_robust_window = 1 in _master.do, then re-run)."
        continue
    }

    local nyr : word count `yrs'
    forvalues i = 1/`nyr' {
        forvalues j = `=`i' + 1'/`nyr' {

            local yi : word `i' of `yrs'
            local yj : word `j' of `yrs'

            use geoid hc_cluster leiden_cluster ///
                using "$temp/final_assignments_temporal_1yr_`yi'_`lvl'.dta", clear
            rename hc_cluster     hc_a
            rename leiden_cluster le_a
            quietly count
            local n_yi = r(N)
            tempfile a
            save `a'

            use geoid hc_cluster leiden_cluster ///
                using "$temp/final_assignments_temporal_1yr_`yj'_`lvl'.dta", clear
            rename hc_cluster     hc_b
            rename leiden_cluster le_b
            quietly count
            local n_yj = r(N)
            merge 1:1 geoid using `a', keep(match) nogenerate
            quietly count
            local n_both = r(N)
            local n_dropped   = `n_yj' - `n_both'
            local n_dropped_a = `n_yi' - `n_both'
            if `n_dropped' > 0 {
                display as text "  [`lvl' `yi'-`yj'] " `n_dropped' ///
                    " of " `n_yj' " `yj'-units absent from `yi' " ///
                    "(compared on the " `n_both' "-unit intersection)"
            }
            if `n_dropped_a' > 0 {
                display as text "  [`lvl' `yi'-`yj'] " `n_dropped_a' ///
                    " of " `n_yi' " `yi'-units absent from `yj' " ///
                    "(compared on the " `n_both' "-unit intersection)"
            }

            foreach m of local methods {
                local va = cond("`m'" == "hc", "hc_a", "le_a")
                local vb = cond("`m'" == "hc", "hc_b", "le_b")

                compare_partitions, a(`va') b(`vb')

                preserve
                clear
                quietly set obs 1
                gen str6   geo_level = "`lvl'"
                gen str8   method    = "`m'"
                gen int    year_a    = `yi'
                gen int    year_b    = `yj'
                gen byte   adjacent  = (`yj' - `yi' == 1)
                gen long   n         = `r(N)'
                gen long   n_dropped = `n_dropped'
                gen double NMI       = `r(NMI)'
                gen double ARI       = `r(ARI)'
                gen double AMI       = `r(AMI)'

                if `yfirst' {
                    save `yout', replace
                    local yfirst 0
                }
                else {
                    append using `yout'
                    save `yout', replace
                }
                restore
            }
        }
    }
}

capture confirm file "`yout'"
if _rc {
    display as error "No year-over-year comparisons produced; single-year runs missing."
    display as error "  Produce them via robust_01/robust_02 behind \$do_robust_window in _master.do."
}
else {
    use `yout', clear
    sort geo_level method year_a year_b

    label variable geo_level "Geography (tract or county)"
    label variable method    "Algorithm"
    label variable year_a    "Earlier year"
    label variable year_b    "Later year"
    label variable adjacent  "1 if consecutive-year pair"
    label variable n         "Units compared (intersection of both years)"
    label variable NMI       "Normalized MI (arithmetic; NOT chance-corrected)"
    label variable ARI       "Adjusted Rand index (chance-corrected)"
    label variable AMI       "Adjusted mutual information (chance-corrected)"
    format NMI ARI AMI %7.4f

    save "$tables/temporal_yoy.dta", replace
    display as result "  -> $tables/temporal_yoy.dta"

    display as text _n "{hline 70}"
    display as text "Year-over-year stability (independent single-year partitions)"
    display as text "{hline 70}"
    list geo_level method year_a year_b adjacent NMI ARI, ///
        noobs sepby(geo_level method)

}

display as result _n "assess_03 done."

cap log close
