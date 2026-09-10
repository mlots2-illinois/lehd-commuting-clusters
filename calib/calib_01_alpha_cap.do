clear all
cap log close
local _logdate = subinstr("$S_DATE", " ", "", .)
log using "$log/calib_01_alpha_cap_${run_tag}_`_logdate'.log", replace text

* calib_01_alpha_cap.do -- sweeps the flow weight alpha and the distance cap over a sample of divisions and writes the dispersion table and figure sources that justify the chosen (alpha, cap) cell.
* Called by: _master.do (calib phase, when $do_calib_sweeps; expects $n_sample_divs, $alpha_grid, $alpha_sc_grid, $cap_grid_<lvl>, $disp_cap_rows_<lvl>, $disp_alpha_rows).
* Reads:  $temp/division_lookup.dta; $temp/division_lookup_county.dta; $temp/pflow_<lvl>_5yr_<ys>_<ye>.dta; includes _inc/sample_divisions.do, _inc/calib_params.do, _inc/build_div_dissim.do
* Writes: $temp/calibration/<lvl>/calib_distance.dta; $temp/calibration/<lvl>/calib_flow.dta; $temp/calibration/<lvl>/calib_sensitivity.dta; $tables/calib_dispersion_<lvl>.csv; $tables/calib_dispersion_grid.csv; $tables/calib_dispersion.tex; via _py/calib_dispersion_grid.py, partition_quality_alpha.py, calib_alpha_view.py, calib_flow_coverage.py and calib_cap_view.py: $tables/partition_quality_alpha_calib_<lvl>.csv, $tables/calib_flow_coverage_<lvl>.csv, $figures/calib_alpha_view_<lvl>.png, $figures/calib_cap_view_<lvl>.png
* Notes:  The chosen alpha and cap come from _inc/calib_params.do, not from this sweep; the sweep documents the dispersion of D_ij around that choice. We keep each unit once, in its core division, before sampling so that overlap units are not counted twice.

include "$program/_inc/sample_divisions.do"
include "$program/_inc/calib_params.do"
include "$program/_inc/build_div_dissim.do"

local n_sample_divs $n_sample_divs
local alpha_grid    "$alpha_grid"
if "`alpha_grid'" == "" | "$cap_grid_tract" == "" | "$cap_grid_county" == "" ///
    | "`n_sample_divs'" == "" {
    display as error "calib_01: \$alpha_grid / \$cap_grid_tract / \$cap_grid_county / \$n_sample_divs not set; run via _master.do."
    exit 198
}

local alpha_sc_grid "$alpha_sc_grid"
if "`alpha_sc_grid'" == "" {
    display as error "calib_01: \$alpha_sc_grid not set; run via _master.do."
    exit 198
}
local alpha_sc_csv = subinstr(trim(itrim("`alpha_sc_grid'")), " ", ",", .)
local pyalpha "$program/_py/partition_quality_alpha.py"
capture confirm file "`pyalpha'"
if _rc {
    display as error "calib_01: missing `pyalpha'; cannot run the alpha selector sweep (Phase D)."
    exit 601
}

capture mkdir "$temp/calibration"

local disp_pieces ""

foreach lvl of global levels {

    display as text _n "{hline 70}"
    display as text "===== LEVEL: `lvl' ====="
    display as text "{hline 70}"

    local outdir "$temp/calibration/`lvl'"
    capture mkdir "`outdir'"

    local pflow_file "$temp/pflow_`lvl'_5yr_${yr_start}_${yr_end}.dta"

    local cap_grid "${cap_grid_`lvl'}"

    local chosen_a  = ``lvl'_alpha'
    local chosen_c  = ``lvl'_cap_km'

    if "`lvl'" == "tract" {
        local divlookup "$temp/division_lookup.dta"
    }
    else {
        local divlookup "$temp/division_lookup_county.dta"
    }
    use div_id geoid lat lon is_core using "`divlookup'", clear

    gen byte _neg_core = -is_core
    bysort geoid (_neg_core div_id): keep if _n == 1
    drop is_core _neg_core

    sample_divisions, n(`n_sample_divs')
    local sample_divs "`r(sample_divs)'"
    local ndivs       `r(ndivs)'
    display as text "Divisions in lookup (`lvl'): `ndivs'"
    display as text "Sample divisions (`: word count `sample_divs'' of `ndivs'): `sample_divs'"

    tempfile divdata
    save `divdata'

    display as text _n "{hline 60}"
    display as text "PHASE A (`lvl'): distance distribution"
    display as text "{hline 60}"

    tempfile pairs_all
    build_div_dissim,                                ///
        divlookup("`divdata'") pflow("`pflow_file'") ///
        restrict_divs("`sample_divs'")               ///
        pairsout("`pairs_all'")

    use `pairs_all', clear

    display as text _n "Pooled distance distribution (`lvl', sampled divisions)"
    display as text "  Total pairs: " _N
    tabstat d_km, s(n mean sd p10 p25 p50 p75 p90 p95 p99 max) ///
        columns(statistics) format(%9.1f)

    display as text _n "Distance percentiles by division:"
    tabstat d_km, by(div_id) ///
        s(n mean p25 p50 p75 p90 p99) ///
        columns(statistics) format(%9.1f) nototal

    display as text _n "Fraction of pairs beyond candidate caps:"
    foreach cap of local cap_grid {
        quietly count if d_km > `cap'
        display as text "  d_km > `cap' km: " %6.2f r(N) / _N * 100 "%"
    }

    preserve
    keep div_id geo_i geo_j d_km
    save "`outdir'/calib_distance.dta", replace
    restore

    display as text _n "{hline 60}"
    display as text "PHASE B (`lvl'): flow distribution (5yr pool)"
    display as text "{hline 60}"

    quietly count if P_ij > 0
    display as text "  Sampled pairs:    " %10.0fc _N
    display as text "  Pairs with flow:  " %10.0fc r(N) ///
        "  (" %5.2f r(N)/_N*100 "%)"

    display as text _n "P_ij (flow > 0 only):"
    tabstat P_ij if P_ij > 0, ///
        s(n mean sd p10 p25 p50 p75 p90 p95 p99 max) ///
        columns(statistics) format(%9.6f)

    gen byte dist_bucket = .
    replace dist_bucket = 1 if d_km <=  10
    replace dist_bucket = 2 if d_km >  10 & d_km <=  25
    replace dist_bucket = 3 if d_km >  25 & d_km <=  50
    replace dist_bucket = 4 if d_km >  50 & d_km <= 100
    replace dist_bucket = 5 if d_km > 100 & d_km <= 200
    replace dist_bucket = 6 if d_km > 200 & d_km <= 500
    replace dist_bucket = 7 if d_km > 500
    label define dbkt 1 "0-10" 2 "10-25" 3 "25-50" 4 "50-100" ///
        5 "100-200" 6 "200-500" 7 "500+", replace
    label values dist_bucket dbkt

    display as text _n "Mean P_ij by distance bucket:"
    tabstat P_ij, by(dist_bucket) s(n mean p50 p90 max) ///
        columns(statistics) format(%9.6f) nototal

    drop dist_bucket
    compress
    save "`outdir'/calib_flow.dta", replace

    display as text _n "{hline 60}"
    display as text "PHASE C (`lvl'): dispersion of D_ij over the (alpha x cap) grid"
    display as text "{hline 60}"

    local grid_py "$program/_py/calib_dispersion_grid.py"
    capture confirm file "`grid_py'"
    if _rc {
        display as error "calib_01: missing `grid_py'; cannot build dispersion grid."
        exit 601
    }
    local alpha_csv = subinstr(trim(itrim("`alpha_grid'")), " ", ",", .)
    local cap_csv   = subinstr(trim(itrim("`cap_grid'")), " ", ",", .)
    local sens_csv  "`outdir'/calib_sensitivity.csv"
    capture erase "`sens_csv'"
    display as text "Calling: $py `grid_py' (`: word count `alpha_grid'' alphas x `: word count `cap_grid'' caps)"
    shell "$py" "`grid_py'" "`outdir'/calib_flow.dta" "`alpha_csv'" "`cap_csv'" "`sens_csv'"
    capture confirm file "`sens_csv'"
    if _rc {
        display as error "calib_01: Phase C helper produced no output (`sens_csv'); check $py env and the log."
        exit 601
    }

    import delimited using "`sens_csv'", varnames(1) clear
    sort alpha cap_km div_id
    compress
    save "`outdir'/calib_sensitivity.dta", replace
    capture erase "`sens_csv'"

    preserve
    collapse (mean) d_mean d_sd iqr d_frac_lo d_frac_hi, by(alpha cap_km)
    sort alpha cap_km
    gen str6 geo_level = "`lvl'"
    order geo_level alpha cap_km d_mean d_sd iqr d_frac_lo d_frac_hi

    label variable d_mean    "Mean D_ij (avg over sampled divisions)"
    label variable d_sd      "SD of D_ij (within-division, avg over divisions)"
    label variable iqr       "IQR of D_ij (avg over divisions)"
    label variable d_frac_lo "Share of pairs with D_ij < 0.05"
    label variable d_frac_hi "Share of pairs with D_ij > 0.95"

    format d_mean d_sd iqr d_frac_lo d_frac_hi %7.4f
    display as text _n "Dispersion of D_ij by (alpha, cap_km) (`lvl'):"
    list alpha cap_km d_mean d_sd iqr d_frac_lo d_frac_hi, ///
        separator(`: word count `cap_grid'') noobs

    export delimited using "$tables/calib_dispersion_`lvl'.csv", replace

    tempfile _disp_`lvl'
    save "`_disp_`lvl''"
    local disp_pieces `"`disp_pieces' "`_disp_`lvl''""'
    restore

    preserve
    use "`outdir'/calib_sensitivity.dta", clear
    collapse (mean) d_mean d_sd iqr d_frac_lo d_frac_hi, by(alpha cap_km)
    quietly summarize d_mean if abs(alpha - `chosen_a') < 1e-6 & cap_km == `chosen_c', meanonly
    local cell_mean = r(mean)
    quietly summarize d_sd if abs(alpha - `chosen_a') < 1e-6 & cap_km == `chosen_c', meanonly
    local cell_sd = r(mean)
    display as result _n "Chosen cell (`lvl'): alpha=`chosen_a' cap_km=`chosen_c'" ///
        "  ->  mean D_ij=" %5.3f `cell_mean' "  SD D_ij=" %5.3f `cell_sd'
    restore

    display as text _n "{hline 60}"
    display as text "PHASE D (`lvl'): self-containment selector vs alpha"
    display as text "{hline 60}"

    local sc_k_nat = ``lvl'_hc_k_national'
    local sc_res   = ``lvl'_leiden_res'
    local samp_csv = subinstr(trim(itrim("`sample_divs'")), " ", ",", .)
    local sc_out   "$tables/partition_quality_alpha_calib_`lvl'.csv"

    display as text "alpha grid: `alpha_sc_grid'  (cap=`chosen_c' k_nat=`sc_k_nat')"
    capture erase "`sc_out'"
    shell PQ_LIGHT=1 PQ_SAMPLE_DIVS="`samp_csv'" "$py" "`pyalpha'" "`lvl'" ///
        "`divlookup'" "`pflow_file'" "`chosen_c'" "`sc_k_nat'" "`sc_res'" ///
        "$seed" "`alpha_sc_csv'" "`sc_out'"
    capture confirm file "`sc_out'"
    if _rc display as error ///
        "  calib_01 Phase D: no selector output for `lvl'; alpha view will fall back."
    else    display as result ///
        "  -> `sc_out' (self-containment vs alpha, sampled divisions)"

    display as result _n "calib_01 (`lvl') diagnostics done. Outputs in `outdir'/."
}

display as text _n "{hline 70}"
display as text "PAPER ARTIFACTS: dispersion table + figure"
display as text "{hline 70}"

clear
append using `disp_pieces'
sort geo_level cap_km alpha
export delimited using "$tables/calib_dispersion_grid.csv", replace
display as result "  -> $tables/calib_dispersion_grid.csv (figure source)"

local n_col = 5
local diags "d_mean d_sd iqr d_frac_hi"
capture file close _tex
file open _tex using "$tables/calib_dispersion.tex", write replace
file write _tex "% Auto-generated by calib_01_alpha_cap.do; do not edit by hand." _n
file write _tex "\begin{tabular}{@{}lrrrr@{}}" _n
file write _tex "\toprule" _n
file write _tex "& Mean \(D_{ij}\) & SD & IQR & Share \(>0.95\) \\" _n
file write _tex "\midrule" _n

local _blank = 0
foreach lvl of global levels {
    local chosen_a  = ``lvl'_alpha'
    local chosen_c  = ``lvl'_cap_km'

    foreach _s in cap alpha {
        if "`_s'" == "cap" local _rows "${disp_cap_rows_`lvl'}"
        else               local _rows "${disp_alpha_rows}"
        if "`_s'" == "cap" local _pin = `chosen_c'
        else               local _pin = `chosen_a'
        local _has = 0
        foreach r of local _rows {
            if abs(`r' - `_pin') < 1e-9 local _has = 1
        }
        if !`_has' local _rows "`_rows' `_pin'"
        numlist "`_rows'", sort
        local `_s'_rows "`r(numlist)'"
    }

    if "`lvl'" != "`: word 1 of $levels'" file write _tex "\addlinespace" _n
    file write _tex "\multicolumn{`n_col'}{@{}l}{\textit{`=strproper("`lvl'")'}} \\[0.15em]" _n

    local af : display %4.2f `chosen_a'
    file write _tex "\multicolumn{`n_col'}{@{}l}{\quad distance cap " ///
        "\(\bar d\) (km), at \(\alpha=`=trim("`af'")'\)} \\" _n
    foreach cap of local cap_rows {
        local capf : display %9.0fc `cap'
        local bold = abs(`cap' - `chosen_c') < 1e-9
        file write _tex "\quad `=trim("`capf'")'"
        foreach st of local diags {
            quietly summarize `st' if geo_level == "`lvl'" & ///
                abs(alpha - `chosen_a') < 1e-6 & cap_km == `cap', meanonly
            if r(N) == 0 {
                file write _tex " & ---"
                local _blank = `_blank' + 1
            }
            else {
                local v : display %5.3f r(mean)
                local v = trim("`v'")
                if `bold' file write _tex " & \textbf{`v'}"
                else      file write _tex " & `v'"
            }
        }
        file write _tex " \\" _n
    }

    file write _tex "\addlinespace[0.3em]" _n
    local capf : display %9.0fc `chosen_c'
    local capm = subinstr(trim("`capf'"), ",", "{,}", .)
    file write _tex "\multicolumn{`n_col'}{@{}l}{\quad flow weight " ///
        "\(\alpha\), at \(\bar d = `capm'\) km} \\" _n
    foreach a of local alpha_rows {
        local af : display %4.2f `a'
        local bold = abs(`a' - `chosen_a') < 1e-9
        file write _tex "\quad `=trim("`af'")'"
        foreach st of local diags {
            quietly summarize `st' if geo_level == "`lvl'" & ///
                abs(alpha - `a') < 1e-6 & cap_km == `chosen_c', meanonly
            if r(N) == 0 {
                file write _tex " & ---"
                local _blank = `_blank' + 1
            }
            else {
                local v : display %5.3f r(mean)
                local v = trim("`v'")
                if `bold' file write _tex " & \textbf{`v'}"
                else      file write _tex " & `v'"
            }
        }
        file write _tex " \\" _n
    }
}
file write _tex "\bottomrule" _n
file write _tex "\end{tabular}" _n
file close _tex
if `_blank' > 0 display as error ///
    "  calib_01: `_blank' dispersion cells not in the sweep; check " ///
    "$disp_cap_rows_* / $disp_alpha_rows against the sweep grids."
display as result "  -> calib_dispersion.tex (LaTeX fragment for tab:calib_dispersion)"

local av_py  "$program/_py/calib_alpha_view.py"
capture confirm file "`av_py'"
if !_rc {
    foreach lvl of global levels {
        local _ca ``lvl'_alpha'
        local _cc ``lvl'_cap_km'
        local sc_csv "$tables/partition_quality_alpha_calib_`lvl'.csv"
        capture confirm file "`sc_csv'"
        if _rc local sc_csv "$tables/partition_quality_alpha_${run_tag}.csv"
        capture confirm file "`sc_csv'"
        if !_rc shell "$py" "`av_py'" "$tables" "$figures" "`lvl'" `_cc' `_ca' "`sc_csv'"
        else    shell "$py" "`av_py'" "$tables" "$figures" "`lvl'" `_cc' `_ca'
        capture confirm file "$figures/calib_alpha_view_`lvl'.png"
        if !_rc display as result ///
            "  -> calib_alpha_view_`lvl'.png (α-on-x: selector vs guard)"
        else    display as error ///
            "  calib_alpha_view.py produced no PNG for `lvl'."
    }
}
else display as error "  Missing `av_py'; skipping alpha companion view."

local fc_py "$program/_py/calib_flow_coverage.py"
local cv_py "$program/_py/calib_cap_view.py"
capture confirm file "`fc_py'"
local _fc_ok = (_rc == 0)
capture confirm file "`cv_py'"
local _cv_ok = (_rc == 0)
if `_fc_ok' & `_cv_ok' {
    foreach lvl of global levels {
        local _ca ``lvl'_alpha'
        local _cc ``lvl'_cap_km'
        local cap_grid "${cap_grid_`lvl'}"
        local cap_csv = subinstr(trim(itrim("`cap_grid'")), " ", ",", .)
        local cov_csv "$tables/calib_flow_coverage_`lvl'.csv"
        local flow_dta "$temp/calibration/`lvl'/calib_flow.dta"

        capture confirm file "`flow_dta'"
        if _rc {
            display as error "  Missing `flow_dta'; skipping cap view for `lvl'."
            continue
        }
        shell "$py" "`fc_py'" "`flow_dta'" "`cap_csv'" "`cov_csv'"

        shell "$py" "`cv_py'" "$tables" "$figures" "`lvl'" `_cc' `_ca'
        capture confirm file "$figures/calib_cap_view_`lvl'.png"
        if !_rc display as result ///
            "  -> calib_cap_view_`lvl'.png (cap selector vs guard)"
        else    display as error ///
            "  calib_cap_view.py produced no PNG for `lvl'."
    }
}
else display as error "  Missing cap-view helpers; skipping cap selector view."

cap log close
