clear all
cap log close
set more off
local _logdate = subinstr("$S_DATE", " ", "", .)
log using "$log/data_02_create_files_${run_tag}_`_logdate'.log", replace text

* data_02_create_files.do -- reads the raw LODES origin-destination (OD) files for every state and year and writes one national tract-level and one county-level OD file per year.
* Called by: _master.do (data phase; once per job type in $job_types, with $job_type and $pool_tag set).
* Reads:  $temp/state_lookup.dta; $raw_lodes/<st>_od_main_<jt>_<year>.csv[.gz]; $raw_lodes/<st>_od_aux_<jt>_<year>.csv[.gz]
* Writes: $temp/1yr/od_s000_<lvl>_1yr_<year><pool_tag>.dta; $log/data_02_coverage_<jt>_<run_tag>_<date><pool_tag>.csv
* Notes:  We append each state's aux file (workers who live out of state) to its main file so that interstate commutes appear in both the home and the work state. Block geocodes are left-padded to 15 characters before the tract (11) and county (5) prefixes are cut. A state short of the full year range is reported in the coverage file, not dropped.

local rawdir "$raw_lodes"
capture mkdir "$temp/1yr"

if "$levels" == "" global levels "tract county"
local do_tract  = strpos("$levels", "tract")  > 0
local do_county = strpos("$levels", "county") > 0

capture program drop import_od_file
program define import_od_file, rclass
    syntax , file(string)
    capture confirm file "`file'.gz"
    if !_rc {
        tempfile tmp
        shell gunzip -c "`file'.gz" > "`tmp'"
        import delimited using "`tmp'", varnames(1) colrange(1:3) stringcols(1 2) clear
        if _N == 0 {
            di as error "import_od_file: `file'.gz extracted to an empty file; corrupt download?"
            exit 459
        }
        rename *, lower
        return scalar found = 1
        exit
    }
    capture confirm file "`file'"
    if !_rc {
        import delimited using "`file'", varnames(1) colrange(1:3) stringcols(1 2) clear
        rename *, lower
        return scalar found = 1
        exit
    }
    return scalar found = 0
end

capture confirm file "$temp/state_lookup.dta"
if _rc {
    di as error "Missing $temp/state_lookup.dta; run data_01 first."
    exit 601
}
use "$temp/state_lookup.dta", clear
levelsof st_abbrev, local(states) clean

foreach st of local states {
    local cov_`st' ""
}

forvalues y = $yr_start/$yr_end {

    display as text _n "{hline 70}"
    display as text "===== data_02 YEAR: `y' ====="
    display as text "{hline 70}"

    local tract_pieces  ""
    local county_pieces ""
    local have_year 0

    foreach st of local states {

        local stem "`rawdir'/`st'_od_main_${job_type}_`y'.csv"
        import_od_file, file("`stem'")
        if !r(found) {
            display as text "  [`y' `st'] no main file; skipping state."
            continue
        }
        keep w_geocode h_geocode s000
        tempfile _main
        quietly save "`_main'"

        local auxstem "`rawdir'/`st'_od_aux_${job_type}_`y'.csv"
        import_od_file, file("`auxstem'")
        if r(found) {
            keep w_geocode h_geocode s000
            append using "`_main'"
        }
        else {
            use "`_main'", clear
            display as text "  [`y' `st'] no aux file; main only."
        }

        capture confirm numeric variable s000
        if _rc {
            di as error "data_02: s000 imported non-numeric for `st' `y'; corrupt LODES file?"
            exit 459
        }
        capture assert !mi(s000)
        if _rc {
            di as error "data_02: missing s000 values for `st' `y'; corrupt LODES file?"
            exit 459
        }
        foreach v in w_geocode h_geocode {
            replace `v' = substr("000000000000000", 1, 15 - length(`v')) + `v' ///
                if length(`v') < 15
        }
        gen str11 w_tract  = substr(w_geocode, 1, 11)
        gen str11 h_tract  = substr(h_geocode, 1, 11)
        gen str5  w_county = substr(w_geocode, 1,  5)
        gen str5  h_county = substr(h_geocode, 1,  5)

        if `do_tract' {
            preserve
            gcollapse (sum) s000, by(w_tract h_tract)
            tempfile _tr_`st'
            quietly save "`_tr_`st''"
            local tract_pieces `"`tract_pieces' "`_tr_`st''""'
            restore
        }
        if `do_county' {
            gcollapse (sum) s000, by(w_county h_county)
            tempfile _co_`st'
            quietly save "`_co_`st''"
            local county_pieces `"`county_pieces' "`_co_`st''""'
        }

        local have_year 1
        local cov_`st' "`cov_`st'' `y'"
        display as text "  [`y' `st'] imported."
    }

    if !`have_year' {
        display as text "  No raw files found for `y'; no 1yr output written."
        continue
    }

    if `do_tract' {
        clear
        append using `tract_pieces'
        gcollapse (sum) s000, by(w_tract h_tract)
        gen int year = `y'
        compress
        save "$temp/1yr/od_s000_tract_1yr_`y'${pool_tag}.dta", replace
        display as result "  Saved od_s000_tract_1yr_`y'${pool_tag}.dta (" _N " pairs)"
    }

    if `do_county' {
        clear
        append using `county_pieces'
        gcollapse (sum) s000, by(w_county h_county)
        gen int year = `y'
        compress
        save "$temp/1yr/od_s000_county_1yr_`y'${pool_tag}.dta", replace
        display as result "  Saved od_s000_county_1yr_`y'${pool_tag}.dta (" _N " pairs)"
    }
}

local exp_years ""
local n_exp = 0
forvalues y = $yr_start/$yr_end {
    local exp_years "`exp_years' `y'"
    local n_exp = `n_exp' + 1
}

display as text _n "{hline 70}"
display as text "===== data_02 STATE-YEAR COVERAGE (${job_type}) ====="
display as text "  Expected `n_exp' year(s) per state:`exp_years'"
display as text "{hline 70}"

tempname covpost
tempfile covtmp
postfile `covpost' str2 state int year byte found using "`covtmp'", replace

local n_short = 0
foreach st of local states {
    local nfound : word count `cov_`st''

    local missing ""
    foreach y of local exp_years {
        local hit = strpos(" `cov_`st'' ", " `y' ") > 0
        post `covpost' ("`st'") (`y') (`hit')
        if !`hit' local missing "`missing' `y'"
    }

    if `nfound' == `n_exp' {
        display as text "  " upper("`st'") ": `nfound'/`n_exp' years; complete."
    }
    else {
        display as error ///
            "  WARNING: " upper("`st'") " has `nfound'/`n_exp' years (missing`missing')"
        local n_short = `n_short' + 1
    }
}
postclose `covpost'

if `n_short' == 0 {
    display as result _n "  All `: word count `states'' states have full `n_exp'-year coverage."
}
else {
    display as error _n "  `n_short' state(s) short of full `n_exp'-year coverage; see warnings above."
}

preserve
use "`covtmp'", clear
export delimited using ///
    "$log/data_02_coverage_${job_type}_${run_tag}_`_logdate'${pool_tag}.csv", replace
restore

display as result _n "data_02 done."

cap log close
