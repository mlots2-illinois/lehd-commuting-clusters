clear all
cap log close
local _logdate = subinstr("$S_DATE", " ", "", .)
log using "$log/data_03_pool_multi_${run_tag}_`_logdate'.log", replace text

* data_03_pool_multi.do -- pools the single-year OD files into 5-year and 3-year windows ending in $yr_end, averaging S000 over the years actually read.
* Called by: _master.do (data phase; once per job type in $job_types, with $job_type and $pool_tag set).
* Reads:  $temp/1yr/od_s000_<lvl>_1yr_<year><pool_tag>.dta
* Writes: $temp/<p>yr/od_s000_<lvl>_<p>yr_<ystart>_<yend><pool_tag>.dta (p = 5, 3)
* Notes:  A pair absent from a year counts as zero in that year. We divide the summed flow by eff_p, the number of window years in which the work state was read, rather than by the window length; in this way a state with a missing LODES year keeps a mean on the same scale as the others.

local pool_windows 5 3
local y_end = $yr_end

foreach p in `pool_windows' {
    capture mkdir "$temp/`p'yr"
}

foreach geo of global levels {

    local wgeo "w_`geo'"
    local hgeo "h_`geo'"

    foreach p in `pool_windows' {

        local y_start = `y_end' - `p' + 1

        display as text _n "Pooling `geo'-level OD: `p'yr (`y_start' - `y_end')"

        use "$temp/1yr/od_s000_`geo'_1yr_`y_start'${pool_tag}.dta", clear
        capture drop year
        rename s000 s000_y0

        local lag = 0
        forvalues y = `=`y_start' + 1'/`y_end' {
            local lag = `lag' + 1
            merge 1:1 `wgeo' `hgeo' using ///
                "$temp/1yr/od_s000_`geo'_1yr_`y'${pool_tag}.dta", ///
                nogen keepusing(s000)
            rename s000 s000_y`lag'
        }

        forvalues k = 0/`=`p' - 1' {
            replace s000_y`k' = 0 if mi(s000_y`k')
        }

        tempfile _pres
        local _first 1
        forvalues y = `y_start'/`y_end' {
            preserve
            use `wgeo' using "$temp/1yr/od_s000_`geo'_1yr_`y'${pool_tag}.dta", clear
            gen str2 wst = substr(`wgeo', 1, 2)
            gcontract wst
            keep wst
            if `_first' {
                quietly save `_pres', replace
                local _first 0
            }
            else {
                append using `_pres'
                quietly save `_pres', replace
            }
            restore
        }
        preserve
        use `_pres', clear
        gcontract wst
        rename _freq eff_p
        label variable eff_p "Years of this window in which the work state was read"
        quietly count if eff_p < `p'
        if r(N) {
            display as error ///
                "  data_03: `r(N)' work state(s) short of the `p'yr window; pooling each over its own years:"
            list wst eff_p if eff_p < `p', noobs clean
        }
        quietly save `_pres', replace
        restore

        gen str2 wst = substr(`wgeo', 1, 2)
        merge m:1 wst using `_pres', keep(match master) nogenerate
        replace eff_p = `p' if mi(eff_p)
        assert eff_p > 0 & eff_p <= `p'

        egen double s000 = rowtotal(s000_y0-s000_y`=`p' - 1')
        replace s000 = s000 / eff_p

        keep `wgeo' `hgeo' s000
        gen int year_start = `y_start'
        gen int year_end   = `y_end'

        label variable s000       "Mean S000 across `p'yr window"
        label variable year_start "First year in pool"
        label variable year_end   "Last year in pool"

        compress
        save "$temp/`p'yr/od_s000_`geo'_`p'yr_`y_start'_`y_end'${pool_tag}.dta", replace
        display as result "  Saved: `p'yr/od_s000_`geo'_`p'yr_`y_start'_`y_end'${pool_tag}.dta (" _N " pairs)"
    }
}

cap log close
