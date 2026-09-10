clear all
cap log close
local _logdate = subinstr("$S_DATE", " ", "", .)
log using "$log/data_04_cbsa_crosswalk_${run_tag}_`_logdate'.log", replace text

* data_04_cbsa_crosswalk.do -- builds county-to-CBSA and tract-to-CBSA crosswalks from the 2023 OMB delineation file for the CBSA comparison in assess_05.
* Called by: _master.do (data phase, once; again at the end of the build phase, after build_02 has written centroids_tract.dta).
* Reads:  $raw/cbsa/list1_2023.xlsx (sheet "List 1"); $temp/centroids_tract.dta
* Writes: $proj/03_data/02_processed/cbsa/county_to_cbsa.dta; $proj/03_data/02_processed/cbsa/tract_to_cbsa.dta
* Notes:  CBSA = core-based statistical area. A tract inherits the CBSA of the county in the first five digits of its GEOID, so the tract file needs centroids_tract.dta; the first call in the data phase writes the county file alone. We drop the states in $noncontig_fips to match the flow universe.

local raw_xlsx "$raw/cbsa/list1_2023.xlsx"
local out_dir  "$proj/03_data/02_processed/cbsa"

capture confirm file "`raw_xlsx'"
if _rc {
    display as error "data_04: missing OMB delineation `raw_xlsx'; cannot build CBSA crosswalk."
    exit 601
}
capture mkdir "$proj/03_data/02_processed"
capture mkdir "`out_dir'"

import excel using "`raw_xlsx'", sheet("List 1") cellrange(A4) allstring clear

keep A D E J K
rename A cbsa
rename D cbsa_title
rename E cbsa_type
rename J stfips
rename K cofips

drop if missing(cbsa) | missing(stfips) | missing(cofips)
keep if regexm(cbsa, "^[0-9]+$")

gen str5 geoid = string(real(stfips), "%02.0f") + string(real(cofips), "%03.0f")
replace cbsa = trim(cbsa)
replace cbsa_title = trim(cbsa_title)
replace cbsa_type  = trim(cbsa_type)

gen str2 _st2 = substr(geoid, 1, 2)
foreach f of global noncontig_fips {
    drop if _st2 == "`f'"
}
drop _st2 stfips cofips

duplicates drop
isid geoid

order geoid cbsa cbsa_title cbsa_type
label variable geoid      "County FIPS (5-digit)"
label variable cbsa       "CBSA code (5-digit)"
label variable cbsa_title "CBSA title"
label variable cbsa_type  "Metropolitan / Micropolitan"
compress
save "`out_dir'/county_to_cbsa.dta", replace

quietly count
display as result "data_04: county_to_cbsa.dta built; `r(N)' counties in a CBSA."
quietly levelsof cbsa, local(_cb)
display as result "         distinct CBSAs: `: word count `_cb''"

tempfile county_xw
save "`county_xw'"

capture confirm file "$temp/centroids_tract.dta"
if _rc {
    display as text "data_04: no $temp/centroids_tract.dta; skipping tract_to_cbsa " ///
        "(rebuild after build_02 when running tract)."
}
else {
    use geoid using "$temp/centroids_tract.dta", clear
    duplicates drop
    gen str5 _county = substr(geoid, 1, 5)
    rename geoid _tract
    rename _county geoid
    merge m:1 geoid using "`county_xw'", keep(match) ///
        keepusing(cbsa cbsa_title cbsa_type) nogenerate
    drop geoid
    rename _tract geoid

    order geoid cbsa cbsa_title cbsa_type
    label variable geoid      "Tract GEOID (11-digit)"
    label variable cbsa       "CBSA code (5-digit)"
    label variable cbsa_title "CBSA title"
    label variable cbsa_type  "Metropolitan / Micropolitan"
    isid geoid
    compress
    save "`out_dir'/tract_to_cbsa.dta", replace

    quietly count
    display as result "data_04: tract_to_cbsa.dta built; `r(N)' tracts in a CBSA."
}

display as result _n "data_04 done. CBSA crosswalk(s) in `out_dir'/."
cap log close
