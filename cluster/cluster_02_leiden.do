clear all
cap log close
local _logdate = subinstr("$S_DATE", " ", "", .)
log using "$log/cluster_02_leiden_${run_tag}_`_logdate'.log", replace text

* cluster_02_leiden.do -- runs Leiden with the constant Potts model (CPM) on the division graphs from cluster_01 at the chosen gamma and at eleven multiples of it.
* Called by: _master.do (cluster phase, $run_tag = baseline) and _inc/run_cluster.do (robust phase; $run_tag, $div_tag, $period_pflow_suffix and $hc_input set by the calling robust_* script).
* Reads:  $temp/chosen_params_<lvl>.dta; $hc_input/<lvl>_<run_tag>/dissim_div<div>.csv; includes _inc/ensure_county_div_lookup.do
* Writes: via _py/leiden_run.py: $temp/leiden/<lvl>/assignments_<run_tag>.csv; $temp/leiden_assignments_<run_tag>_<lvl>.dta
* Notes:  The multiples (0.10 to 3.00 times the chosen gamma) feed the gamma sensitivity in the assess and gen phases; cluster_03 keeps only the resolution nearest the chosen one. Every run uses $seed.

include "$program/_inc/ensure_county_div_lookup.do"

local pyscript "$program/_py/leiden_run.py"

foreach lvl of global levels {

    display as text _n "{hline 70}"
    display as text "===== cluster_02 LEVEL: `lvl' ====="
    display as text "{hline 70}"

    local params     "$temp/chosen_params_`lvl'.dta"
    local dissim_dir "$hc_input/`lvl'_${run_tag}"

    local outdir     "$temp/leiden/`lvl'"
    capture mkdir "$temp/leiden"
    capture mkdir "`outdir'"
    local result_csv "`outdir'/assignments_${run_tag}.csv"

    capture confirm file "`dissim_dir'/dissim_div00001.csv"
    if _rc {
        local _anycsv : dir "`dissim_dir'" files "dissim_div*.csv"
        if `"`_anycsv'"' == "" {
            display as error "cluster_02: no dissim CSVs in `dissim_dir'; run cluster_01 first."
            exit 601
        }
    }

    use "`params'", clear
    local res_chosen = leiden_resolution[1]
    local res_grid ""
    foreach mult in 0.10 0.20 0.35 0.50 0.70 0.85 1.00 1.20 1.50 2.00 2.50 3.00 {
        local res_grid `res_grid' `=`res_chosen' * `mult''
    }
    local res_csv : subinstr local res_grid " " ",", all
    if substr("`res_csv'", 1, 1) == "," local res_csv = substr("`res_csv'", 2, .)
    display as text "Leiden resolutions (`lvl'): `res_grid'"

    display as text _n "Calling: $py `pyscript' `dissim_dir' ... (cpm)"
    capture erase "`result_csv'"
    shell "$py" "`pyscript'" "`dissim_dir'" "`result_csv'" "`res_csv'" "$seed" "cpm"
    capture confirm file "`result_csv'"
    if _rc {
        display as error "cluster_02: Python step produced no output (`result_csv'); check $py env and the log."
        exit 601
    }

    import delimited using "`result_csv'", varnames(1) ///
        stringcols(2) clear
    confirm string variable geoid
    assert strlen(geoid) == cond("`lvl'" == "tract", 11, 5)

    label variable div_id          "Division ID"
    label variable geoid           "Geographic FIPS (`lvl')"
    label variable resolution      "Leiden resolution parameter"
    label variable leiden_cluster  "Cluster ID within (div_id, resolution)"

    compress
    save "$temp/leiden_assignments_${run_tag}_`lvl'.dta", replace

    display as result _n ///
        "cluster_02 (`lvl') done. Saved leiden_assignments_${run_tag}_`lvl'.dta: " _N " rows."
}

display as result _n "cluster_02 done for $levels."

cap log close
