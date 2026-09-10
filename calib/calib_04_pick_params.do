clear all
cap log close
local _logdate = subinstr("$S_DATE", " ", "", .)
log using "$log/calib_04_pick_params_${run_tag}_`_logdate'.log", replace text

* calib_04_pick_params.do -- persists the hand-chosen alpha, distance cap, HC k and Leiden gamma from _inc/calib_params.do as one-row files the cluster scripts read.
* Called by: _master.do (calib phase; once, after the sweeps).
* Reads:  includes _inc/calib_params.do; $temp/chosen_params_<lvl>.dta (read back for display)
* Writes: $temp/chosen_params_<lvl>.dta
* Notes:  The cluster and robust scripts take parameters only from these files. robust_06 and robust_08 overwrite them per specification and restore from _chosen_params_baseline_<lvl>.dta; we erase any such backup here because the file just written is the baseline again.

include "$program/_inc/calib_params.do"

foreach lvl of global levels {

    local leiden_res = ``lvl'_leiden_res'

    clear
    set obs 1
    gen str6   geo_level         = "`lvl'"
    gen double alpha             = ``lvl'_alpha'
    gen double cap_km            = ``lvl'_cap_km'
    gen long   hc_k_national     = ``lvl'_hc_k_national'
    gen double leiden_resolution = `leiden_res'

    label variable geo_level         "Geography (tract or county)"
    label variable alpha             "Flow vs. distance weight in D_ij"
    label variable cap_km            "Distance cap (km)"
    label variable hc_k_national     "National target HC cluster count"
    label variable leiden_resolution "Leiden resolution parameter"

    save "$temp/chosen_params_`lvl'.dta", replace

    * A backup left by a crashed robust_06 or robust_08 run is stale once chosen_params is rewritten from calib_params.do.
    capture erase "$temp/_chosen_params_baseline_`lvl'.dta"

    display as result _n "Chosen parameters (`lvl'):"
    list, noobs
}

* These formatted locals are not used below; gen/tab_chosen_params.do reads the .dta files and formats the table itself.
foreach lvl in tract county {
    local s = cond("`lvl'" == "tract", "t", "c")
    local a_`s' "--"
    local c_`s' "--"
    local k_`s' "--"
    local g_`s' "--"
    capture confirm file "$temp/chosen_params_`lvl'.dta"
    if !_rc & strpos("$levels", "`lvl'") {
        use "$temp/chosen_params_`lvl'.dta", clear
        local a_`s' : display %4.2f alpha[1]
        local c_`s' : display %5.0f cap_km[1]
        local k_`s' : display %9.0fc hc_k_national[1]
        local g_`s' : display %5.2f leiden_resolution[1]
        local k_`s' = subinstr(trim("`k_`s''"), ",", "{,}", .)
    }
    else if !strpos("$levels", "`lvl'") {
        display as text "  [`lvl'] not in \$levels; chosen_params column rendered as --."
    }
}

display as result "  -> chosen_params_<lvl>.dta persisted " ///
    "(table fragment is emitted by gen/tab_chosen_params.do in the GENERATE phase)"

cap log close
