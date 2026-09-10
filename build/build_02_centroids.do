clear all
cap log close
local _logdate = subinstr("$S_DATE", " ", "", .)
log using "$log/build_02_centroids_${run_tag}_`_logdate'.log", replace text

* build_02_centroids.do -- builds the tract and county centroid files from the Census 2020 tract centers of population.
* Called by: _master.do (build phase; once).
* Reads:  $raw/tract_populationcenters_2020.csv
* Writes: $temp/centroids_tract.dta; $temp/centroids_county.dta
* Notes:  County centroids are population-weighted means of their tract centroids, not the Census county centers of population. Tracts in $noncontig_fips and tracts with missing coordinates are dropped.

import delimited using "$raw/tract_populationcenters_2020.csv", ///
    delimiters(",") varnames(1) stringcols(1 2 3) clear

rename (statefp countyfp tractce population latitude longitude) ///
    (st_fips co_fips tr_fips pop lat lon)

destring st_fips, generate(_st_num)
foreach _fp of global noncontig_fips {
    drop if _st_num == real("`_fp'")
}
drop _st_num

gen str11 geoid = st_fips + co_fips + tr_fips
gen str5  county_geoid = st_fips + co_fips

drop if mi(lat) | mi(lon)

label variable geoid        "11-digit tract FIPS"
label variable county_geoid "5-digit county FIPS"
label variable st_fips      "State FIPS"
label variable lat          "Population center latitude"
label variable lon          "Population center longitude"
label variable pop          "Tract population (2020)"

preserve
keep geoid st_fips lat lon pop
order geoid st_fips lat lon pop
sort geoid
compress
save "$temp/centroids_tract.dta", replace
display as result "  centroids_tract.dta:  " _N " tracts"
restore

collapse (mean) lat lon [aweight = pop], by(county_geoid st_fips)
rename county_geoid geoid

label variable geoid   "5-digit county FIPS"
label variable st_fips "State FIPS"
label variable lat     "Pop-weighted county centroid latitude"
label variable lon     "Pop-weighted county centroid longitude"

order geoid st_fips lat lon
sort geoid
compress
save "$temp/centroids_county.dta", replace
display as result "  centroids_county.dta: " _N " counties"

cap log close
