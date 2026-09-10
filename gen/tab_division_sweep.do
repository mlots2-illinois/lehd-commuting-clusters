clear all

* tab_division_sweep.do -- write the division-geometry table comparing each robust_07_divisions.do rerun with the tract baseline, labelled by bands, window W and step S.
* Called by: _master.do (generate phase, do_generate switch).
* Reads:  includes _inc/booktabs.do and _inc/compare_to_baseline.do; $temp/final_assignments_baseline_tract.dta; $temp/final_assignments_div_<div_tag>_tract.dta
* Writes: $tables/division_sweep.tex; $tables/division_sweep.csv
* Notes:  The row list is $div_specs_tract from _inc/div_specs.do, so geometries that were built but not rerun in robust_07 are skipped with a message.

include "$program/_inc/booktabs.do"
include "$program/_inc/compare_to_baseline.do"

capture confirm file "$temp/final_assignments_baseline_tract.dta"
if _rc {
    display as text "tab_division_sweep: no tract baseline; skipping."
    exit 0
}
if `"$div_specs_tract"' == "" {
    display as text "tab_division_sweep: \$div_specs_tract not defined (run via _master.do); skipping."
    exit 0
}

capture program drop _div_label
program define _div_label, rclass
    args tag
    local lab "`tag'"
    if "`tag'" == "baseline"  local lab "Baseline"
    if "`tag'" == "Wnarrow"   local lab "Narrower window"
    if "`tag'" == "Wmid"      local lab "Mid-width window"
    if "`tag'" == "Wwide"     local lab "Wider window"
    if "`tag'" == "Wxwide"    local lab "Extra-wide window"
    if "`tag'" == "ovxsmall"  local lab "Much smaller overlap"
    if "`tag'" == "ovsmall"   local lab "Smaller overlap"
    if "`tag'" == "ovlarge"   local lab "Larger overlap"
    if "`tag'" == "ovxlarge"  local lab "Much larger overlap"
    if "`tag'" == "bands15"   local lab "Coarsest bands"
    if "`tag'" == "bands20"   local lab "Coarser bands"
    if "`tag'" == "bands45"   local lab "Finer bands"
    if "`tag'" == "bands60"   local lab "Finest bands"
    if "`tag'" == "smallboth" local lab "Smaller window and step"
    if "`tag'" == "bigboth"   local lab "Larger window and step"
    return local lab "`lab'"
end

tempname pf
tempfile res
postfile `pf' str12 tag str24 setting int bands int W int S ///
    double NMIhc double ARIhc double NMIle double ARIle using "`res'", replace

foreach s of global div_specs_tract {
    gettoken tag   rest : s
    gettoken bands rest : rest
    gettoken W     S    : rest
    _div_label `tag'
    local setting `"`r(lab)'"'

    if "`tag'" == "baseline" {
        post `pf' ("baseline") (`"`setting'"') (`bands') (`W') (`S') (1) (1) (1) (1)
        continue
    }
    capture confirm file "$temp/final_assignments_div_`tag'_tract.dta"
    if _rc {
        display as text "  no div_`tag' rerun; skipping that row."
        continue
    }
    compare_to_baseline, lvl(tract) alt("$temp/final_assignments_div_`tag'_tract.dta")
    local nh = r(NMI_hc)
    local ah = r(ARI_hc)
    local nl = r(NMI_le)
    local al = r(ARI_le)
    post `pf' ("`tag'") (`"`setting'"') (`bands') (`W') (`S') (`nh') (`ah') (`nl') (`al')
}
postclose `pf'

use "`res'", clear
quietly count
if r(N) == 0 {
    display as text "tab_division_sweep: nothing to write; skipping."
    exit 0
}
export delimited tag setting bands W S NMIhc ARIhc NMIle ARIle ///
    using "$tables/division_sweep.csv", replace

bt_open, handle(_t) path("$tables/division_sweep.tex") colspec(l ccc cc cc) ///
    script(gen/tab_division_sweep.do)
file write _t "& & & & \multicolumn{2}{c}{\textbf{Hierarchical}} & \multicolumn{2}{c}{\textbf{Leiden}} \\" _n
file write _t "\cmidrule(lr){5-6} \cmidrule(lr){7-8}" _n
file write _t "\textbf{Setting} & \textbf{Bands} & \textbf{\(W\)} & \textbf{\(S\)} & NMI & ARI & NMI & ARI \\" _n
file write _t "\midrule" _n
forvalues i = 1/`=_N' {
    tex_escape, s(`=setting[`i']')
    local set `"`r(s)'"'
    if tag[`i'] == "baseline" local set "`set'\,$^\dagger$"
    bt_num, value(`=W[`i']') fmt(%5.0fc)
    local wv `"`r(s)'"'
    bt_num, value(`=S[`i']') fmt(%5.0fc)
    local sv `"`r(s)'"'
    bt_num, value(`=NMIhc[`i']') fmt(%5.3f)
    local nh `"`r(s)'"'
    bt_num, value(`=ARIhc[`i']') fmt(%5.3f)
    local ah `"`r(s)'"'
    bt_num, value(`=NMIle[`i']') fmt(%5.3f)
    local nl `"`r(s)'"'
    bt_num, value(`=ARIle[`i']') fmt(%5.3f)
    local al `"`r(s)'"'
    bt_row, handle(_t) cells(`" "`set'" "`=bands[`i']'" "`wv'" "`sv'" "`nh'" "`ah'" "`nl'" "`al'" "')
    if `i' == 1 file write _t "\addlinespace" _n
}
bt_close, handle(_t)

display as result "  -> division_sweep.tex + division_sweep.csv"
