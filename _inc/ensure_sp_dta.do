* ensure_sp_dta.do -- defines ensure_sp_dta, which converts the 2020 tract, county or state shapefile to Stata sp format once and returns r(ok).
* Included by: gen/fig_map_national.do, fig_map_regions.do, fig_map_states.do, fig_flow_map.do, fig_walk_il.do, fig_walk_tristate.do, fig_walk_omaha.do, fig_walk_chicago.do.
* Expects: $shapefiles, $temp
* Reads:  $shapefiles/us_<level>_2020/US_<level>_2020.shp
* Writes: $temp/<tracts|counties|states>_sp.dta and the paired _sp_shp.dta (spshape2dta)
* Notes:  a missing shapefile returns r(ok) = 0 rather than an error, so the map scripts skip cleanly on a machine without the shapefiles.

capture program drop ensure_sp_dta
program define ensure_sp_dta, rclass

    syntax , level(string)

    if "`level'" == "tract" {
        local sp_basename "tracts_sp"
        local shp_file    "$shapefiles/us_tract_2020/US_tract_2020.shp"
    }
    else if "`level'" == "state" {
        local sp_basename "states_sp"
        local shp_file    "$shapefiles/us_state_2020/US_state_2020.shp"
    }
    else {
        local sp_basename "counties_sp"
        local shp_file    "$shapefiles/us_county_2020/US_county_2020.shp"
    }
    local sp_dta "$temp/`sp_basename'.dta"

    capture confirm file "`sp_dta'"
    if !_rc {
        return scalar ok = 1
        exit 0
    }

    capture confirm file "`shp_file'"
    if _rc {
        di as text "ensure_sp_dta: shapefile not found: `shp_file'"
        di as text "  Cannot build `sp_dta'; place the 2020 `level' shapefile"
        di as text "  under 03_data/00_shapefiles/ and rerun.  Skipping."
        return scalar ok = 0
        exit 0
    }

    di as text "Converting `level' shapefile (one-time) -> `sp_dta' ..."
    local saved_pwd : pwd
    quietly cd "$temp"
    spshape2dta "`shp_file'", saving("`sp_basename'") replace
    quietly cd "`saved_pwd'"
    return scalar ok = 1
end
