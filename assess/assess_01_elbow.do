clear all
cap log close
local _logdate = subinstr("$S_DATE", " ", "", .)
log using "$log/assess_01_elbow_${run_tag}_`_logdate'.log", replace text

* assess_01_elbow.do -- compute the mean within-cluster dissimilarity along an HC k grid and along the Leiden resolution grid, for the elbow figures, and pick the knee of each curve.
* Called by: _master.do (assess phase, do_assess switch).
* Reads:  $hc_input/<lvl>_<run_tag>/dissim_div<div>.csv; $temp/division_lookup[_county].dta; $temp/leiden_assignments_<run_tag>_<lvl>.dta; shell "$py" _py/hc_run.py
* Writes: $temp/elbow_<run_tag>_<lvl>.dta; $temp/elbow_pick_<run_tag>_<lvl>.dta; $temp/_elbow_hc_assign_<k>_<lvl>.csv (erased after import)
* Notes:  The HC grid here is coarser and wider than the calibration grid in calib_03_hc_k.do because the figure needs the whole curve, not the neighbourhood of the chosen k. The Leiden side reuses the sweep assignments already stored by cluster_02_leiden.do, so no Leiden run happens here. The knee is the point farthest from the chord between the curve's endpoints after both axes are scaled to [0, 1].

local k_grid_tract  1000 3000 5000 7500 10000 12500 15000 17500 20000 22500 25000 27500 30000 32500 35000 37500 40000 45000 50000 60000 75000 90000 105000 120000
local k_grid_county 50 100 200 300 400 500 600 700 800 900 1000 1100 1200 1300 1400 1500 1750 2000 2500 3000

local pyscript "$program/_py/hc_run.py"

foreach lvl of global levels {

    display as text _n "{hline 70}"
    display as text "===== assess_01 LEVEL: `lvl' ====="
    display as text "{hline 70}"

    local hc_dir "$hc_input/`lvl'_${run_tag}"
    local _dfiles ""
    capture local _dfiles : dir "`hc_dir'" files "dissim_div*.csv"
    if `:word count `_dfiles'' == 0 {
        display as error "No dissim CSVs in `hc_dir'; skipping `lvl'."
        display as error "  Re-run cluster_01 (run_tag=${run_tag}) to regenerate the cached"
        display as error "  per-division dissimilarity CSVs, then re-run assess_01."
        continue
    }

    if "`lvl'" == "tract" {
        local k_grid    `k_grid_tract'
        local divlookup "$temp/division_lookup.dta"
    }
    else {
        local k_grid    `k_grid_county'
        local divlookup "$temp/division_lookup_county.dta"
    }

    use div_id geoid is_core using "`divlookup'", clear
    gen byte _neg_core = -is_core
    bysort geoid (_neg_core div_id): keep if _n == 1
    drop is_core _neg_core
    quietly count
    local total_n = r(N)

    bysort div_id: gen long n_div_count = _N
    bysort div_id: keep if _n == 1
    keep div_id n_div_count
    tempfile div_sizes
    save `div_sizes'

    capture frame drop _acc
    frame create _acc str8 method double param long div_id long n_clusters ///
        double wcsd_sum double n_within_pairs

    foreach k_nat of local k_grid {
        display as text _n "  --- HC sweep: national k = `k_nat' (running hc_run.py) ---"

        local avg_size = `total_n' / `k_nat'

        use `div_sizes', clear
        gen long k_div = round(n_div_count / `avg_size')
        replace k_div = max(2, k_div)
        replace k_div = min(k_div, n_div_count)

        local k_map_str ""
        forvalues r = 1/`=_N' {
            local d  = div_id[`r']
            local kd = k_div[`r']
            if "`k_map_str'" == "" {
                local k_map_str "`d':`kd'"
            }
            else {
                local k_map_str "`k_map_str',`d':`kd'"
            }
        }

        local out_csv "$temp/_elbow_hc_assign_`k_nat'_`lvl'.csv"
        capture erase "`out_csv'"
        shell "$py" "`pyscript'" "`hc_dir'" "`out_csv'" --k-map "`k_map_str'"
        capture confirm file "`out_csv'"
        if _rc {
            display as error "assess_01: hc_run.py produced no output for k=`k_nat' (`lvl'); check $py env."
            cap log close
            exit 601
        }

        import delimited using "`out_csv'", varnames(1) stringcols(2) clear
        rename hc_cluster cluster
        keep div_id geoid cluster
        tempfile hc_assign_`k_nat'
        save `hc_assign_`k_nat'', replace
        erase "`out_csv'"
    }

    quietly use `div_sizes', clear
    quietly levelsof div_id, local(all_divs)

    foreach d of local all_divs {
        local fname = string(`d', "$divfmt")
        local dissim "`hc_dir'/dissim_div`fname'.csv"
        capture confirm file "`dissim'"
        if _rc continue

        display as text "  --- HC WCSD: division `d' (import dissim once, sweep k) ---"
        quietly import delimited using "`dissim'", varnames(1) ///
            stringcols(1 2) clear
        capture rename d D
        tempfile dissim_d
        quietly save `dissim_d', replace

        foreach k_nat of local k_grid {

            quietly use if div_id == `d' using "`hc_assign_`k_nat''", clear
            quietly count
            if r(N) == 0 continue
            drop div_id

            preserve
            contract cluster
            local nclus = _N
            restore

            rename cluster cluster_i
            rename geoid   geo_i
            tempfile clu_i
            quietly save `clu_i', replace
            rename geo_i     geo_j
            rename cluster_i cluster_j
            tempfile clu_j
            quietly save `clu_j', replace

            quietly use `dissim_d', clear
            quietly merge m:1 geo_i using `clu_i', keep(match) nogenerate
            quietly merge m:1 geo_j using `clu_j', keep(match) nogenerate

            quietly keep if cluster_i == cluster_j
            quietly summarize D
            local wsum = cond(r(N) > 0, r(sum), 0)
            local wnp  = r(N)

            frame post _acc ("HC") (`k_nat') (`d') (`nclus') (`wsum') (`wnp')
        }
    }

    capture frame drop _le
    frame create _le
    frame _le: use "$temp/leiden_assignments_${run_tag}_`lvl'.dta", clear

    frame _le: egen long _res_id = group(resolution)
    frame _le: quietly levelsof _res_id, local(res_ids)
    foreach rid of local res_ids {
        frame _le: quietly summarize resolution if _res_id == `rid', meanonly
        local resval_`rid' = r(mean)
    }
    frame _le: quietly levelsof div_id, local(divs)

    foreach d of local divs {
        local fname = string(`d', "$divfmt")
        local dissim "`hc_dir'/dissim_div`fname'.csv"
        capture confirm file "`dissim'"
        if _rc continue

        quietly import delimited using "`dissim'", varnames(1) ///
            stringcols(1 2) clear
        capture rename d D
        tempfile dissim_d
        quietly save `dissim_d'

        foreach rid of local res_ids {
            quietly use `dissim_d', clear

            frame _le {
                preserve
                quietly keep if div_id == `d' & _res_id == `rid'
                keep geoid leiden_cluster
                tempfile clus
                quietly save `clus'
                restore
            }

            rename geo_i geoid
            quietly merge m:1 geoid using `clus', keep(match) nogenerate
            rename geoid geo_i
            rename leiden_cluster cluster_i

            rename geo_j geoid
            quietly merge m:1 geoid using `clus', keep(match) nogenerate
            rename geoid geo_j
            rename leiden_cluster cluster_j

            quietly keep if cluster_i == cluster_j
            quietly summarize D
            local wsum = cond(r(N) > 0, r(sum), 0)
            local wnp  = r(N)

            frame _le {
                quietly levelsof leiden_cluster ///
                    if div_id == `d' & _res_id == `rid', local(cl)
            }
            local nclus : word count `cl'

            frame post _acc ("Leiden") (`resval_`rid'') (`d') (`nclus') (`wsum') (`wnp')
        }
    }

    frame change _acc
    gcollapse (sum) n_clusters wcsd_sum n_within_pairs, by(method param)
    sort method param

    gen double wcsd_mean = cond(n_within_pairs > 0, wcsd_sum / n_within_pairs, .)

    label variable method         "Algorithm"
    label variable param          "k_national (HC) or resolution (Leiden)"
    label variable n_clusters     "Total clusters across all divisions"
    label variable wcsd_sum       "Within-cluster SUM of dissimilarities (reference)"
    label variable n_within_pairs "Count of within-cluster pairs"
    label variable wcsd_mean      "MEAN within-cluster dissimilarity (elbow Y-axis)"

    compress
    save "$temp/elbow_${run_tag}_`lvl'.dta", replace
    display as result _n "assess_01 (`lvl') done. " _N " rows."
    list method param n_clusters wcsd_mean, noobs sepby(method)

    capture frame drop _pick
    frame create _pick str8 method double picked_param long picked_n_clusters double knee_score

    quietly levelsof method, local(methods)
    foreach mth of local methods {
        preserve
        quietly keep if method == "`mth'" & !mi(wcsd_mean)
        sort n_clusters
        quietly count
        if r(N) >= 3 {
            quietly summarize n_clusters
            local xmin = r(min)
            local xrng = r(max) - r(min)
            quietly summarize wcsd_mean
            local ymin = r(min)
            local yrng = r(max) - r(min)
            if `xrng' > 0 & `yrng' > 0 {
                gen double _xn = (n_clusters - `xmin') / `xrng'
                gen double _yn = (wcsd_mean  - `ymin') / `yrng'
                local x1 = _xn[1]
                local y1 = _yn[1]
                local x2 = _xn[_N]
                local y2 = _yn[_N]
                local dx = `x2' - `x1'
                local dy = `y2' - `y1'
                local den = sqrt((`dx')^2 + (`dy')^2)
                gen double _dist = abs((`dy')*_xn - (`dx')*_yn + (`x2')*(`y1') - (`y2')*(`x1')) / `den'
                gsort -_dist
                local pp = param[1]
                local pn = n_clusters[1]
                local ks = _dist[1]
                frame post _pick ("`mth'") (`pp') (`pn') (`ks')
                display as result "  [`lvl' `mth'] knee at param=`pp' (n_clusters=`pn', score=" %5.3f `ks' ")"
            }
        }
        restore
    }

    frame change _pick
    label variable method            "Algorithm"
    label variable picked_param       "Selected k_national (HC) / resolution (Leiden) at the knee"
    label variable picked_n_clusters  "Realized total clusters at the knee"
    label variable knee_score         "Max normalized distance-to-chord (relative knee strength)"
    capture save "$temp/elbow_pick_${run_tag}_`lvl'.dta", replace

    frame change default
    frame drop _acc
    frame drop _le
    capture frame drop _pick
}

display as result _n "assess_01 done for tract and county."

cap log close
