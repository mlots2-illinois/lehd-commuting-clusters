clear all
cap log close
local _logdate = subinstr("$S_DATE", " ", "", .)
log using "$log/data_01_state_lookup_${run_tag}_`_logdate'.log", replace text

* data_01_state_lookup.do -- builds the lookup of two-digit state FIPS codes and lowercase abbreviations that data_02 loops over when reading LODES files.
* Called by: _master.do (data phase; once).
* Writes: $temp/state_lookup.dta
* Notes:  The list holds the 48 contiguous states and DC. We omit AK (02) and HI (15) here so that no LODES file for either is ever read; $noncontig_fips repeats the exclusion downstream for the territories.

input str2 st_fips str2 st_abbrev
    "01" "al"
    "04" "az"
    "05" "ar"
    "06" "ca"
    "08" "co"
    "09" "ct"
    "10" "de"
    "11" "dc"
    "12" "fl"
    "13" "ga"
    "16" "id"
    "17" "il"
    "18" "in"
    "19" "ia"
    "20" "ks"
    "21" "ky"
    "22" "la"
    "23" "me"
    "24" "md"
    "25" "ma"
    "26" "mi"
    "27" "mn"
    "28" "ms"
    "29" "mo"
    "30" "mt"
    "31" "ne"
    "32" "nv"
    "33" "nh"
    "34" "nj"
    "35" "nm"
    "36" "ny"
    "37" "nc"
    "38" "nd"
    "39" "oh"
    "40" "ok"
    "41" "or"
    "42" "pa"
    "44" "ri"
    "45" "sc"
    "46" "sd"
    "47" "tn"
    "48" "tx"
    "49" "ut"
    "50" "vt"
    "51" "va"
    "53" "wa"
    "54" "wv"
    "55" "wi"
    "56" "wy"
end

label variable st_fips   "2-digit state FIPS"
label variable st_abbrev "Lowercase state abbreviation"

save "$temp/state_lookup.dta", replace

cap log close
