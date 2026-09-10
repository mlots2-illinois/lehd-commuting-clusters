clear all
cap log close
local _logdate = subinstr("$S_DATE", " ", "", .)
log using "$log/calib_03_hc_k_${run_tag}_`_logdate'.log", replace text

* calib_03_hc_k.do -- sweeps the national HC cluster target k over a sample of divisions and writes the size and spread summaries and the validity-vs-k source from which k is chosen.
* Called by: _master.do (calib phase, when $do_calib_sweeps; expects $n_sample_divs, $k_grid_<lvl>).
* Reads:  $temp/division_lookup.dta; $temp/division_lookup_county.dta; $temp/pflow_<lvl>_<period_pflow_suffix>.dta; $temp/centroids_<lvl>.dta; includes _inc/calib_params.do, _inc/sample_divisions.do, _inc/build_div_dissim.do
* Writes: $temp/calibration/<lvl>/hc_calib_input/dissim_div<div>.csv; $temp/calibration/<lvl>/calib_hc_k.dta; $temp/calibration/<lvl>/calib_hc_k_summary.dta; $tables/calib_hc_k_summary_<lvl>.csv; via _py/hc_run.py, partition_quality_alpha.py and calib_hc_k_view.py: $tables/partition_quality_k_<lvl>.csv, $figures/calib_hc_k_view_<lvl>.png
* Notes:  HC = hierarchical clustering. Each national k is allocated to a division in proportion to its unit count, floored at 2 and capped at the division size, the same rule cluster_01 applies. All per-division cuts go to hc_run.py in one call so each dendrogram is built once and cut many times. The bounding-box diagonal is the within-cluster spread guard.

local n_sample_divs $n_sample_divs
local k_grid_tract  "$k_grid_tract"
local k_grid_county "$k_grid_county"
if "`k_grid_tract'" == "" | "`k_grid_county'" == "" | "`n_sample_divs'" == "" {
    display as error "calib_03: \$k_grid_tract / \$k_grid_county / \$n_sample_divs not set; run via _master.do."
    exit 198
}

include "$program/_inc/calib_params.do"
include "$program/_inc/sample_divisions.do"
include "$program/_inc/build_div_dissim.do"

local pyscript "$program/_py/hc_run.py"

capture mkdir "$temp/calibration"

foreach lvl of global levels {

    display as text _n "{hline 70}"
    display as text "===== calib_03 LEVEL: `lvl' ====="
    display as text "{hline 70}"

    local outdir "$temp/calibration/`lvl'"
    capture mkdir "`outdir'"

    if "`lvl'" == "tract" {
        local divlookup "$temp/division_lookup.dta"
        local k_grid    `k_grid_tract'
    }
    else {
        local divlookup "$temp/division_lookup_county.dta"
        local k_grid    `k_grid_county'
    }
    local alpha  ``lvl'_alpha'
    local cap_km ``lvl'_cap_km'

    local asrc "_inc/calib_params.do"
    local csrc "_inc/calib_params.do"
    local pflow_file "$temp/pflow_`lvl'_${period_pflow_suffix}.dta"
    display as text "Using alpha=`alpha' (`asrc')  cap_km=`cap_km' (`csrc')"

    use div_id geoid lat lon is_core using "`divlookup'", clear
    gen byte _neg_core = -is_core
    bysort geoid (_neg_core div_id): keep if _n == 1
    drop is_core _neg_core

    quietly count
    local total_n = r(N)
    quietly levelsof div_id, local(all_divs)
    local ndivs : word count `all_divs'
    display as text "Total `lvl' units: `total_n'    Divisions: `ndivs'"

    sample_divisions, n(`n_sample_divs')
    local sample_divs "`r(sample_divs)'"
    display as text "Sample divisions (`: word count `sample_divs'' of `ndivs'): `sample_divs'"

    tempfile divdata
    save `divdata'

    local hc_dir "`outdir'/hc_calib_input"
    capture mkdir "`hc_dir'"

    build_div_dissim,                                ///
        divlookup("`divdata'") pflow("`pflow_file'") ///
        alpha(`alpha') cap_km(`cap_km')              ///
        outdir("`hc_dir'") divfmt("%05.0f")          ///
        restrict_divs("`sample_divs'")               ///
        countsframe(_c03_counts)

    tempfile sample_sizes
    frame _c03_counts {
        rename n_units n_div
        sort div_id
        save "`sample_sizes'"
    }
    frame drop _c03_counts

    local kalloc_pieces ""
    foreach k_nat of local k_grid {
        use `sample_sizes', clear
        local avg_size = `total_n' / `k_nat'
        gen long k_target_national = `k_nat'
        gen long k_div = round(n_div / `avg_size')
        replace k_div = max(2, k_div)
        replace k_div = min(k_div, n_div)
        keep k_target_national div_id k_div
        tempfile _ka_`k_nat'
        quietly save "`_ka_`k_nat''"
        local kalloc_pieces `"`kalloc_pieces' "`_ka_`k_nat''""'
    }
    clear
    append using `kalloc_pieces'
    sort k_target_national div_id
    tempfile kalloc
    save `kalloc'

    keep div_id k_div
    duplicates drop
    sort div_id k_div

    local k_list_str ""
    quietly levelsof div_id, local(_kl_divs)
    foreach d of local _kl_divs {
        local _ks ""
        quietly levelsof k_div if div_id == `d', local(_kvals) clean
        foreach kv of local _kvals {
            if "`_ks'" == "" {
                local _ks "`kv'"
            }
            else {
                local _ks "`_ks'|`kv'"
            }
        }
        if "`k_list_str'" == "" {
            local k_list_str "`d':`_ks'"
        }
        else {
            local k_list_str "`k_list_str',`d':`_ks'"
        }
    }

    local cache_csv "`outdir'/_hc_k_assign_all.csv"
    capture erase "`cache_csv'"
    shell "$py" "`pyscript'" "`hc_dir'" "`cache_csv'" --k-list "`k_list_str'"
    capture confirm file "`cache_csv'"
    if _rc {
        display as error "calib_03: Python step produced no output (`cache_csv'); check $py env and the log."
        exit 601
    }

    import delimited using "`cache_csv'", varnames(1) stringcols(2) clear
    confirm string variable geoid
    assert strlen(geoid) == cond("`lvl'" == "tract", 11, 5)
    gen long _k_div = round(cut_value)
    drop cut_value
    rename _k_div k_div

    merge m:1 geoid using "$temp/centroids_`lvl'.dta", ///
        keepusing(lat lon) keep(match master) nogenerate
    tempfile hc_cache
    save `hc_cache'

    local hck_pieces ""

    foreach k_nat of local k_grid {
        display as text _n "  --- HC sweep: national k = `k_nat' ---"

        local avg_size = `total_n' / `k_nat'
        display as text "    avg cluster size target = " %6.2f `avg_size' " units"

        use `kalloc' if k_target_national == `k_nat', clear
        keep div_id k_div
        tempfile _sel
        save `_sel'

        use `hc_cache', clear
        merge m:1 div_id k_div using `_sel', keep(match) nogenerate
        gen double k_requested        = k_div
        gen long   k_target_national  = `k_nat'

        bysort div_id hc_cluster: gen long _cluster_n = _N
        bysort div_id hc_cluster: egen double _lat_min = min(lat)
        bysort div_id hc_cluster: egen double _lat_max = max(lat)
        bysort div_id hc_cluster: egen double _lon_min = min(lon)
        bysort div_id hc_cluster: egen double _lon_max = max(lon)
        local torad = _pi / 180
        gen double _bbox_km = 2 * 6371 * asin(sqrt( ///
            sin((_lat_max - _lat_min)*`torad'/2)^2 + ///
            cos(_lat_min*`torad')*cos(_lat_max*`torad') ///
            *sin((_lon_max - _lon_min)*`torad'/2)^2))

        bysort div_id hc_cluster: keep if _n == 1
        keep k_target_national div_id hc_cluster k_requested _cluster_n _bbox_km
        rename _cluster_n cluster_n
        rename _bbox_km   bbox_km

        tempfile _hck_`k_nat'
        quietly save "`_hck_`k_nat''"
        local hck_pieces `"`hck_pieces' "`_hck_`k_nat''""'
    }

    capture erase "`cache_csv'"

    clear
    append using `hck_pieces'
    sort k_target_national div_id hc_cluster
    label variable k_target_national "National k target"
    label variable div_id            "Division ID"
    label variable hc_cluster        "Cluster ID within division"
    label variable k_requested       "k requested for this division"
    label variable cluster_n         "Number of units in cluster"
    label variable bbox_km           "Bounding-box diagonal (km)"
    compress
    save "`outdir'/calib_hc_k.dta", replace

    preserve
    gcollapse                                                      ///
        (count) n_clusters_obs = hc_cluster                         ///
        (mean)  mean_size      = cluster_n  mean_bbox = bbox_km     ///
        (p50)   p50_size       = cluster_n  p50_bbox  = bbox_km     ///
        (max)   max_size       = cluster_n  max_bbox  = bbox_km     ///
        (min)   min_size       = cluster_n                          ///
        , by(k_target_national)
    format mean_size p50_size %7.2f
    format mean_bbox p50_bbox max_bbox %7.1f
    sort k_target_national
    list k_target_national n_clusters_obs ///
        min_size p50_size mean_size max_size ///
        p50_bbox mean_bbox max_bbox, noobs sepby(k_target_national)
    save "`outdir'/calib_hc_k_summary.dta", replace
    export delimited using "$tables/calib_hc_k_summary_`lvl'.csv", replace
    restore

    local pyalpha "$program/_py/partition_quality_alpha.py"
    local kview   "$program/_py/calib_hc_k_view.py"
    capture confirm file "`pyalpha'"
    local _pa_ok = (_rc == 0)
    capture confirm file "`kview'"
    local _kv_ok = (_rc == 0)
    if `_pa_ok' & `_kv_ok' {
        local kgrid_csv = subinstr(trim(itrim("`k_grid'")), " ", ",", .)
        local samp_csv  = subinstr(trim(itrim("`sample_divs'")), " ", ",", .)
        local alpha_csv "`alpha'"
        local chosen_k  = ``lvl'_hc_k_national'
        local res       = ``lvl'_leiden_res'
        local kq_csv "$tables/partition_quality_k_`lvl'.csv"
        capture erase "`kq_csv'"
        display as text "  validity vs k (k grid: `k_grid'; alpha=`alpha' cap=`cap_km')"
        * PQ_FULL=1 keeps silhouette (Rousseeuw 1987) and modularity switched
        * on across the k grid; without it partition_quality_alpha.py runs in
        * light mode and returns them as missing, leaving the cut with no
        * internal-validity evidence at all.
        shell PQ_KGRID="`kgrid_csv'" PQ_SAMPLE_DIVS="`samp_csv'" PQ_FULL=1 "$py" ///
            "`pyalpha'" "`lvl'" "`divlookup'" "`pflow_file'" "`cap_km'" ///
            "`chosen_k'" "`res'" "$seed" "`alpha_csv'" "`kq_csv'"
        capture confirm file "`kq_csv'"
        if _rc {
            display as error "  no self-containment-vs-k output for `lvl'; check \$py env and log."
        }
        else {
            display as result "  -> `kq_csv' (self-containment, silhouette, modularity vs k)"
            shell "$py" "`kview'" "$tables" "$figures" "`lvl'" `chosen_k'
            capture confirm file "$figures/calib_hc_k_view_`lvl'.png"
            if !_rc display as result ///
                "  -> calib_hc_k_view_`lvl'.png (SC bound vs spread guard)"
            else    display as error ///
                "  calib_hc_k_view.py produced no PNG for `lvl'."
        }
    }
    else display as error "  Missing k-criterion helpers; skipping validity-vs-k view."

    display as result _n "calib_03 (`lvl') done. Outputs in `outdir'/."
    display as result "    Choose hc_k_national from: $tables/calib_hc_k_summary_`lvl'.csv"
    display as result "    + figure calib_hc_k_view_`lvl'.png (self-containment bound vs spread guard)"
}

display as result _n "calib_03 done for tract and county. Inspect calib_hc_k_summary.dta to choose hc_k_national."

cap log close
