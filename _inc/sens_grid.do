* sens_grid.do -- builds the one-at-a-time sensitivity grid: each calibrated parameter (alpha, cap_km, hc_k_national, leiden_resolution) at 0.5x, 0.75x, 1.25x and 1.5x its chosen value, per level.
* Included by: robust/robust_06_params.do, assess/assess_04_sensitivity.do.
* Expects: $temp, $levels; defines the locals sens_tags, sens_params, sens_v_t and sens_v_c in the including scope.
* Reads:  $temp/chosen_params_tract.dta; $temp/chosen_params_county.dta (fallback: _inc/calib_params.do)
* Notes:  alpha is clipped to (0, 1] and k is rounded to an integer of at least 2. A multiple that yields the same values as an earlier one (alpha at 1.5x with alpha = 0.8 clips to 1, as 1.25x does) is dropped, so the tags sens_<abbr>_m<mult x 100> are unique.

local mults "0.50 0.75 1.25 1.50"

include "$program/_inc/calib_params.do"
local a_t   `tract_alpha'
local cap_t `tract_cap_km'
local k_t   `tract_hc_k_national'
local res_t `tract_leiden_res'
local a_c   `county_alpha'
local cap_c `county_cap_km'
local k_c   `county_hc_k_national'
local res_c `county_leiden_res'

capture confirm file "$temp/chosen_params_tract.dta"
if !_rc {
    preserve
    use "$temp/chosen_params_tract.dta", clear
    local a_t   = alpha[1]
    local cap_t = cap_km[1]
    local k_t   = hc_k_national[1]
    local res_t = leiden_resolution[1]
    restore
}
capture confirm file "$temp/chosen_params_county.dta"
if !_rc {
    preserve
    use "$temp/chosen_params_county.dta", clear
    local a_c   = alpha[1]
    local cap_c = cap_km[1]
    local k_c   = hc_k_national[1]
    local res_c = leiden_resolution[1]
    restore
}

local sens_tags   ""
local sens_params ""
local sens_v_t    ""
local sens_v_c    ""

foreach row in "alpha a_t a_c a" "cap_km cap_t cap_c cap" ///
    "hc_k_national k_t k_c k" "leiden_resolution res_t res_c res" {
    gettoken pname row : row
    gettoken bt    row : row
    gettoken bc    row : row
    gettoken pabb  row : row
    local base_t = ``bt''
    local base_c = ``bc''

    local seen ""

    foreach m of local mults {
        local vt = `base_t' * `m'
        local vc = `base_c' * `m'

        if "`pname'" == "alpha" {
            local vt = min(max(`vt', 0.0001), 1)
            local vc = min(max(`vc', 0.0001), 1)
        }
        else if "`pname'" == "hc_k_national" {
            local vt = max(round(`vt'), 2)
            local vc = max(round(`vc'), 2)
        }

        local key ""
        foreach L of global levels {
            local key "`key'|`=cond("`L'"=="tract", "`vt'", "`vc'")'"
        }
        local dup : list key in seen
        if `dup' continue
        local seen `"`seen' `key'"'

        local mtag = string(`m'*100, "%03.0f")
        local sens_tags   `"`sens_tags' sens_`pabb'_m`mtag'"'
        local sens_params `"`sens_params' `pname'"'
        local sens_v_t    `"`sens_v_t' `vt'"'
        local sens_v_c    `"`sens_v_c' `vc'"'
    }
}
