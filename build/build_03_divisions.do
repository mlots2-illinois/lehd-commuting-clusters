clear all
cap log close

local dtag "$div_tag"
local tagsuf = cond("`dtag'" == "baseline", "", "_`dtag'")

local _logdate = subinstr("$S_DATE", " ", "", .)
log using "$log/build_03_divisions_${run_tag}`tagsuf'_`_logdate'.log", replace text

* build_03_divisions.do -- partitions the tract or county centroids into overlapping snake-scan divisions, the windows within which HC and Leiden run.
* Called by: _master.do (build phase; once per level and per spec in $div_specs_<lvl>, with $div_tag, $div_bands, $div_levels and $div_W/$div_S (tract) or $div_W_county/$div_S_county set).
* Reads:  $temp/centroids_<lvl>.dta
* Writes: $temp/division_lookup<tag>.dta; $temp/division_meta<tag>.dta; $temp/division_lookup_county<tag>.dta; $temp/division_meta_county<tag>.dta (<tag> empty for baseline, else _<div_tag>)
* Notes:  HC = hierarchical clustering. Units are sorted into $div_bands latitude bands and read in alternating directions from band to band, so consecutive ranks stay spatially close at the band edges. A division is a window of W consecutive ranks stepped by S; the W - S units shared with the next window form the overlap zone (is_core = 0) through which cluster_04 fuses clusters across divisions.

local n_bands = $div_bands

local levels "$div_levels"
if "`levels'" == "" local levels "tract county"

foreach lvl of local levels {

    display as text _n "{hline 70}"
    display as text "===== LEVEL: `lvl' ====="
    display as text "{hline 70}"

    if "`lvl'" == "tract" {
        local W         = $div_W
        local S         = $div_S
        local centroids "$temp/centroids_tract.dta"
        local lkp_out   "$temp/division_lookup`tagsuf'.dta"
        local meta_out  "$temp/division_meta`tagsuf'.dta"
    }
    else {
        local W         = $div_W_county
        local S         = $div_S_county
        local centroids "$temp/centroids_county.dta"
        local lkp_out   "$temp/division_lookup_county`tagsuf'.dta"
        local meta_out  "$temp/division_meta_county`tagsuf'.dta"
    }

    use geoid st_fips lat lon using "`centroids'", clear
    local N = _N
    display as text "Loaded `N' `lvl' units."

    xtile lat_band = lat, nq(`n_bands')
    gen double lon_sort = cond(mod(lat_band, 2) == 1, -lon, lon)
    gsort -lat_band lon_sort
    gen long spatial_rank = _n
    drop lon_sort

    label variable spatial_rank "Position in snake-scan spatial sort"
    label variable lat_band     "Latitude band (xtile by lat; 1 = southernmost)"

    local overlap = `W' - `S'
    local ndivs = ceil((`N' - `W') / `S') + 1

    display as text ///
        "Window=`W' step=`S' overlap=`overlap' -> `ndivs' divisions for `N' units."

    tempfile units
    sort spatial_rank
    quietly save `units'

    clear
    quietly set obs `ndivs'
    gen long div_id = _n
    gen long lo     = (div_id - 1) * `S' + 1
    gen long hi     = min(lo + `W' - 1, `N')
    gen long n_in   = hi - lo + 1

    expand n_in
    bysort div_id (lo): gen long offset = _n - 1
    gen long spatial_rank = lo + offset
    keep div_id spatial_rank

    merge m:1 spatial_rank using `units', ///
        assert(match using) keep(match) nogenerate

    quietly {
        gen long div_lo = (div_id - 1) * `S' + 1
        gen long div_hi = min(div_lo + `W' - 1, `N')
        gen long dist_from_lo = spatial_rank - div_lo
        gen long dist_from_hi = div_hi - spatial_rank

        gen byte is_core = ///
            (dist_from_lo >= `overlap' | div_id == 1) & ///
            (dist_from_hi >= `overlap' | div_id == `ndivs')

        label define corelbl 0 "Overlap" 1 "Core", replace
        label values is_core corelbl
        drop div_lo div_hi dist_from_lo dist_from_hi
    }

    destring st_fips, generate(statefp)

    order div_id geoid statefp lat lon spatial_rank lat_band is_core
    sort  div_id spatial_rank

    label variable div_id   "Division ID"
    label variable geoid    "Geography FIPS"
    label variable statefp  "State FIPS"
    label variable lat      "Centroid latitude"
    label variable lon      "Centroid longitude"
    label variable is_core  "1 = core of window, 0 = overlap zone"

    compress
    save "`lkp_out'", replace

    local NL = _N
    display as result _n "Saved `lkp_out'; `NL' rows (`N' units x 1-2 divisions)."

    preserve
    collapse (count) n_units  = spatial_rank ///
        (sum)   n_core   = is_core      ///
        (min)   lat_min  = lat          ///
        (max)   lat_max  = lat          ///
        (min)   lon_min  = lon          ///
        (max)   lon_max  = lon          ///
        (min)   rank_min = spatial_rank ///
        (max)   rank_max = spatial_rank ///
        , by(div_id)

    gen int n_overlap = n_units - n_core

    label variable n_units   "Total units in division"
    label variable n_core    "Units in core (non-overlap) zone"
    label variable n_overlap "Units in overlap zone"

    order div_id n_units n_core n_overlap lat_min lat_max lon_min lon_max
    sort div_id

    save "`meta_out'", replace

    quietly summarize n_units
    display as text   _n "{hline 60}"
    display as text   "Snake-scan division summary (`lvl')"
    display as text   "{hline 60}"
    display as result "  Divisions:        " `ndivs'
    display as result "  Units/div min:    " %5.0f r(min)
    display as result "  Units/div max:    " %5.0f r(max)
    display as result "  Units/div mean:   " %5.0f r(mean)
    quietly summarize n_overlap
    display as result "  Overlap mean:     " %5.0f r(mean)
    display as result "  Total rows:       " `NL'
    display as result "  Avg div/unit:     " %4.2f (`NL' / `N')
    restore
}

cap log close
