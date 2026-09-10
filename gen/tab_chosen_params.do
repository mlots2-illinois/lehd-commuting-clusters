clear all

* tab_chosen_params.do -- write the table of calibrated parameters (alpha, dbar, k, gamma) for both geographies, with the cluster counts the baseline run delivered.
* Called by: _master.do (generate phase, do_generate switch).
* Reads:  includes _inc/booktabs.do and _inc/get_chosen.do; $temp/chosen_params_<lvl>.dta; $temp/final_assignments_<run_tag>_<lvl>.dta
* Writes: $tables/chosen_params.tex; $tables/chosen_params.csv
* Notes:  Delivered counts are taken after seam fusion and contiguity repair, so they need not equal the target k.

include "$program/_inc/booktabs.do"
include "$program/_inc/get_chosen.do"

capture confirm file "$temp/chosen_params_tract.dta"
local have_t = !_rc
capture confirm file "$temp/chosen_params_county.dta"
local have_c = !_rc
if !`have_t' & !`have_c' {
    display as text "tab_chosen_params: no chosen_params files (run calib_04); skipping."
    exit 0
}

foreach lvl in tract county {
    get_chosen, lvl(`lvl') param(alpha)
    local a_`lvl' = r(value)
    get_chosen, lvl(`lvl') param(cap_km)
    local c_`lvl' = r(value)
    get_chosen, lvl(`lvl') param(hc_k_national)
    local k_`lvl' = r(value)
    get_chosen, lvl(`lvl') param(leiden_resolution)
    local g_`lvl' = r(value)
}

foreach lvl in tract county {
    local nhc_`lvl' = .
    local nle_`lvl' = .
    capture confirm file "$temp/final_assignments_${run_tag}_`lvl'.dta"
    if !_rc {
        foreach m in hc leiden {
            preserve
            use `m'_cluster using "$temp/final_assignments_${run_tag}_`lvl'.dta", clear
            quietly drop if mi(`m'_cluster)
            gcontract `m'_cluster
            if "`m'" == "hc" local nhc_`lvl' = _N
            else             local nle_`lvl' = _N
            restore
        }
    }
    else display as text ///
        "tab_chosen_params: no final_assignments_${run_tag}_`lvl'.dta; delivered counts left blank."
}

preserve
clear
set obs 6
gen str20 parameter = ""
gen str8  symbol    = ""
gen double tract    = .
gen double county   = .
replace parameter = "flow_distance_weight" in 1
replace symbol = "alpha"  in 1
replace tract  = `a_tract'  in 1
replace county = `a_county' in 1
replace parameter = "distance_cap_km" in 2
replace symbol = "dbar"   in 2
replace tract  = `c_tract'  in 2
replace county = `c_county' in 2
replace parameter = "hc_target_k" in 3
replace symbol = "k"      in 3
replace tract  = `k_tract'  in 3
replace county = `k_county' in 3
replace parameter = "leiden_resolution" in 4
replace symbol = "gamma"  in 4
replace tract  = `g_tract'  in 4
replace county = `g_county' in 4
replace parameter = "clusters_delivered_hc" in 5
replace symbol = "n_HC"   in 5
replace tract  = `nhc_tract'  in 5
replace county = `nhc_county' in 5
replace parameter = "clusters_delivered_leiden" in 6
replace symbol = "n_Leiden" in 6
replace tract  = `nle_tract'  in 6
replace county = `nle_county' in 6
export delimited using "$tables/chosen_params.csv", replace
restore

bt_num, value(`a_tract') fmt(%4.2f)
local at `"`r(s)'"'
bt_num, value(`a_county') fmt(%4.2f)
local ac `"`r(s)'"'
bt_num, value(`c_tract') fmt(%5.0f)
local ct `"`r(s)'"'
bt_num, value(`c_county') fmt(%5.0f)
local cc `"`r(s)'"'
bt_num, value(`k_tract') fmt(%9.0fc)
local kt `"`r(s)'"'
bt_num, value(`k_county') fmt(%9.0fc)
local kc `"`r(s)'"'
bt_num, value(`g_tract') fmt(%5.3f)
local gt `"`r(s)'"'
bt_num, value(`g_county') fmt(%5.3f)
local gc `"`r(s)'"'

bt_open, handle(_t) path("$tables/chosen_params.tex") colspec(l c r r l) ///
    script(gen/tab_chosen_params.do)
file write _t "Parameter & Symbol & Tract & County & Selection rule \\" _n
file write _t "\midrule" _n
bt_row, handle(_t) cells(`" "Flow--distance weight" "\(\alpha\)" "`at'" "`ac'" "dispersion of \(D_{ij}\) (\texttt{calib\_01})" "')
bt_row, handle(_t) cells(`" "Distance cap (km)" "\(\bar{d}\)" "`ct'" "`cc'" "dispersion + commuting-horizon (\texttt{calib\_01})" "')
bt_row, handle(_t) cells(`" "HC target cluster count" "\(k\)" "`kt'" "`kc'" "over-extension share + fragmentation floor (\texttt{calib\_03})" "')
bt_row, handle(_t) cells(`" "Leiden resolution" "\(\gamma\)" "`gt'" "`gc'" "finest reproducible non-degenerate (\texttt{calib\_02})" "')

bt_num, value(`nhc_tract')  fmt(%9.0fc)
local dht `"`r(s)'"'
bt_num, value(`nhc_county') fmt(%9.0fc)
local dhc `"`r(s)'"'
bt_num, value(`nle_tract')  fmt(%9.0fc)
local dlt `"`r(s)'"'
bt_num, value(`nle_county') fmt(%9.0fc)
local dlc `"`r(s)'"'
file write _t "\addlinespace" _n
file write _t "\multicolumn{5}{l}{\textit{Delivered by the baseline run}} \\" _n
bt_row, handle(_t) cells(`" "Clusters, hierarchical" "" "`dht'" "`dhc'" "after seam fusion and contiguity repair" "')
bt_row, handle(_t) cells(`" "Clusters, Leiden" "" "`dlt'" "`dlc'" "after seam fusion and contiguity repair" "')
bt_close, handle(_t)

display as result "  -> chosen_params.tex + chosen_params.csv"
