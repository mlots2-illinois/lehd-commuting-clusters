clear all
cap log close
local _logdate = subinstr("$S_DATE", " ", "", .)
log using "$log/calib_02_resolution_${run_tag}_`_logdate'.log", replace text

* calib_02_resolution.do -- sweeps the Leiden CPM resolution gamma over a sample of divisions and writes the sweep table and line graph from which gamma is chosen by hand.
* Called by: _master.do (calib phase, when $do_calib_sweeps; expects $n_sample_divs, $resolutions_<lvl>, $leiden_nseeds_<lvl>).
* Reads:  $temp/division_lookup.dta; $temp/division_lookup_county.dta; $temp/pflow_<lvl>_5yr_<ys>_<ye>.dta; includes _inc/build_div_dissim.do, _inc/calib_params.do, _inc/sample_divisions.do
* Writes: $temp/calibration/<lvl>/dissim/dissim_div<div>.csv; $temp/calibration/<lvl>/calib_resolution.dta; $temp/calibration/<lvl>/leiden_resolution_sweep.dta; $tables/leiden_resolution_sweep_<lvl>.csv; via _py/leiden_calibrate.py and leiden_resolution_lines.py: $figures/leiden_resolution_lines_<lvl>.png
* Notes:  CPM = constant Potts model. The national cluster count is the sampled sum rescaled by the ratio of all divisions to sampled divisions, so it is an estimate. The picked gamma is then written by hand into _inc/calib_params.do; this script reads it back only to mark it on the figure.

local n_sample_divs $n_sample_divs
if "`n_sample_divs'" == "" {
    display as error "calib_02: \$n_sample_divs not set; run via _master.do."
    exit 198
}

local leiden_ptype        cpm
local resolutions_tract   "$resolutions_tract"
local resolutions_county  "$resolutions_county"
if "`resolutions_tract'" == "" | "`resolutions_county'" == "" {
    display as error "calib_02: \$resolutions_tract / \$resolutions_county not set; run via _master.do."
    exit 198
}

local pyscript      "$program/_py/leiden_calibrate.py"

capture mkdir "$temp/calibration"

include "$program/_inc/build_div_dissim.do"
include "$program/_inc/calib_params.do"
include "$program/_inc/sample_divisions.do"

foreach lvl of global levels {

    display as text _n "{hline 70}"
    display as text "===== LEVEL: `lvl' ====="
    display as text "{hline 70}"

    local outdir "$temp/calibration/`lvl'"
    capture mkdir "`outdir'"

    capture erase "`outdir'/picked_leiden_res.dta"

    local pflow_file "$temp/pflow_`lvl'_5yr_${yr_start}_${yr_end}.dta"

    if "`lvl'" == "tract" {
        local divlookup "$temp/division_lookup.dta"
        local resolutions "`resolutions_tract'"
    }
    else {
        local divlookup "$temp/division_lookup_county.dta"
        local resolutions "`resolutions_county'"
    }

    local nseeds "${leiden_nseeds_`lvl'}"
    if "`nseeds'" == "" local nseeds 5

    use div_id geoid using "`divlookup'", clear

    sample_divisions, n(`n_sample_divs')
    local sample_divs "`r(sample_divs)'"
    local ndivs       `r(ndivs)'
    display as text "Sample divisions for Leiden calibration (`lvl', `: word count `sample_divs'' of `ndivs'): `sample_divs'"

    local alpha  = ``lvl'_alpha'
    local cap_km = ``lvl'_cap_km'
    display as text "Dense graph (`lvl'): alpha=`alpha' cap_km=`cap_km'"

    local dissim_dir "`outdir'/dissim"
    build_div_dissim,                              ///
        divlookup("`divlookup'")                   ///
        pflow("`pflow_file'")                      ///
        alpha(`alpha') cap_km(`cap_km')            ///
        outdir("`dissim_dir'") divfmt("$divfmt")   ///
        restrict_divs("`sample_divs'")

    local out_csv "`outdir'/calib_resolution.csv"

    display as text _n "Calling: $py `pyscript' `dissim_dir' ... (`leiden_ptype', `nseeds' seeds)"
    capture erase "`out_csv'"
    shell LEIDEN_NSEEDS=`nseeds' "$py" "`pyscript'" "`dissim_dir'" "`out_csv'" "`resolutions'" "$seed" "`leiden_ptype'"
    capture confirm file "`out_csv'"
    if _rc {
        display as error "calib_02: Python step produced no output (`out_csv'); check $py env and the log."
        exit 601
    }

    import delimited using "`out_csv'", varnames(1) clear

    label variable div_id     "Division ID"
    label variable resolution "Leiden CPM resolution (gamma)"
    label variable n_clusters "Mean number of communities (across seeds)"
    label variable modularity "Mean graph modularity (reference)"
    capture label variable stability    "Cross-seed stability (mean pairwise NMI)"
    capture label variable largest_frac "Mean largest-community fraction"
    capture label variable singleton_frac "Mean singleton vertex fraction"

    compress
    save "`outdir'/calib_resolution.dta", replace

    display as text _n "n_clusters by (gamma, division) [`lvl']:"
    tabdisp resolution div_id, cellvar(n_clusters)

    display as text _n "Mean across divisions (`lvl'): n_clusters, stability, largest_frac, singleton_frac:"
    preserve
    collapse (mean) n_clusters stability largest_frac singleton_frac, by(resolution)
    format n_clusters %7.1f
    format stability largest_frac singleton_frac %7.4f
    list resolution n_clusters stability largest_frac singleton_frac, noobs
    restore

    use div_id using "`divlookup'", clear
    duplicates drop
    quietly count
    local n_divs_total = r(N)
    local n_divs_sampled : word count `sample_divs'

    use "`outdir'/calib_resolution.dta", clear

    quietly duplicates tag div_id resolution, gen(_dup)
    quietly count if _dup > 0
    if r(N) {
        display as error "  calib_02: `r(N)' duplicated (division, gamma) row(s) " ///
            "in the `lvl' sweep; keeping one of each before summarising."
        quietly duplicates drop div_id resolution, force
    }
    drop _dup

    preserve
    collapse (sum)  n_clusters_sampled = n_clusters     ///
        (mean) mean_stability     = stability       ///
        mean_largest_frac  = largest_frac    ///
        mean_singleton_frac = singleton_frac ///
        mean_modularity    = modularity,     ///
        by(resolution)
    gen double n_clusters_est_national = ///
        n_clusters_sampled * (`n_divs_total' / `n_divs_sampled')

    label variable resolution              "Leiden CPM resolution (gamma)"
    label variable n_clusters_sampled      "Sum of n_clusters across `n_divs_sampled' sampled divs"
    label variable n_clusters_est_national "Estimated national n_clusters (linear rescale)"
    label variable mean_stability          "Mean cross-seed stability across sampled divs"
    label variable mean_largest_frac       "Mean largest-community fraction (coarse-end degeneracy)"
    label variable mean_singleton_frac     "Mean singleton vertex fraction (fine-end degeneracy)"
    label variable mean_modularity         "Mean graph modularity (reference)"

    format n_clusters_est_national %12.0f
    format mean_stability mean_largest_frac mean_singleton_frac mean_modularity %7.4f
    sort resolution

    save "`outdir'/leiden_resolution_sweep.dta", replace
    export delimited using "$tables/leiden_resolution_sweep_`lvl'.csv", replace

    display as text _n "Sweep summary (`lvl'):"
    list resolution n_clusters_sampled n_clusters_est_national ///
        mean_stability mean_largest_frac mean_singleton_frac, noobs

    restore

    local fig_py "$program/_py/leiden_resolution_lines.py"
    capture confirm file "`fig_py'"
    if !_rc {
        local _cr ``lvl'_leiden_res'
        shell "$py" "`fig_py'" "$tables" "$figures" "`lvl'" --picked-gamma `_cr' --compact
        capture confirm file "$figures/leiden_resolution_lines_`lvl'.png"
        if !_rc display as result "    -> leiden_resolution_lines_`lvl'.png (selector vs guard)"
        else     display as error "    leiden_resolution_lines.py ran but no PNG (`lvl')."

    }
    else display as error "    Missing `fig_py'; skipping Leiden line graph (`lvl')."

    display as result _n "calib_02 (`lvl') done. Outputs in `outdir'/."
    display as result "    Sweep table: $tables/leiden_resolution_sweep_`lvl'.csv"
    display as result "    Line graph:  $figures/leiden_resolution_lines_`lvl'.png"
}

display as result _n "calib_02 done for tract and county. Sweep tables + line graphs produced; Set tract_leiden_res / county_leiden_res by hand in _inc/calib_params.do."

cap log close
