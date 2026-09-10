clear all
cap log close
local _logdate = subinstr("$S_DATE", " ", "", .)
log using "$log/fig_map_states_${run_tag}_`_logdate'.log", replace text

* fig_map_states.do -- draw a per-state atlas of the tract partition, one map per state and method, showing clusters that cross the state line.
* Called by: _master.do (generate phase, do_generate switch).
* Reads:  includes _inc/ensure_sp_dta.do and _inc/state_hangover_map.do; via state_hangover_map, $temp/tracts_sp.dta, $temp/states_sp.dta and $temp/final_assignments_<run_tag>_tract.dta
* Writes: via state_hangover_map, $figures/map_state_<abbr>_<method>_<run_tag>.pdf
* Notes:  The list holds the 48 conterminous states plus DC in FIPS order; Alaska, Hawaii and the territories are excluded upstream by $noncontig_fips.

include "$program/_inc/ensure_sp_dta.do"
include "$program/_inc/state_hangover_map.do"

local fips "01 04 05 06 08 09 10 11 12 13 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 44 45 46 47 48 49 50 51 53 54 55 56"
local abbr "AL AZ AR CA CO CT DE DC FL GA ID IL IN IA KS KY LA ME MD MA MI MN MS MO MT NE NV NH NJ NM NY NC ND OH OK OR PA RI SC SD TN TX UT VT VA WA WV WI WY"

local n : word count `fips'
forvalues i = 1/`n' {
    local f : word `i' of `fips'
    local a : word `i' of `abbr'
    display as text _n "{hline 70}"
    display as text "fig_map_states: `a' (FIPS `f')  [`i'/`n']"
    display as text "{hline 70}"
    state_hangover_map, state("`f'") suffix("`a'")
}

display as result _n "fig_map_states done.  Per-state atlas written to $figures/."

cap log close
