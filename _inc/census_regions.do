* census_regions.do -- state FIPS lists for the four Census regions and the nine Census divisions, plus the key lists the callers loop over.
* Included by: gen/fig_map_regions.do (regions); assess/assess_09_maup.do, assess/assess_11_robinson.do, assess/assess_12_random_null.do (divisions).
* Expects: nothing.
* Notes:  Alaska and Hawaii are listed in $div_pacific for completeness; $noncontig_fips removes them before any unit is built, so the Pacific division is CA, OR and WA in practice.

global states_northeast "09 23 25 33 34 36 42 44 50"

global states_midwest   "17 18 19 20 26 27 29 31 38 39 46 55"

global states_south     "01 05 10 11 12 13 21 22 24 28 37 40 45 47 48 51 54"

global states_west      "04 06 08 16 30 32 35 41 49 53 56"

global region_keys "northeast midwest south west"

* Census divisions (nine), one rung of the areal-unit ladder in assess_09, assess_11 and assess_12.
global div_newengland  "09 23 25 33 44 50"
global div_midatlantic "34 36 42"
global div_encentral   "17 18 26 39 55"
global div_wncentral   "19 20 27 29 31 38 46"
global div_satlantic   "10 11 12 13 24 37 45 51 54"
global div_escentral   "01 21 28 47"
global div_wscentral   "05 22 40 48"
global div_mountain    "04 08 16 30 32 35 49 56"
global div_pacific     "02 06 15 41 53"

global division_keys "newengland midatlantic encentral wncentral satlantic escentral wscentral mountain pacific"
