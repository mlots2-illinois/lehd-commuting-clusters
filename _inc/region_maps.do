* region_maps.do -- defines region_maps, which draws tract-level HC and Leiden maps for a set of focal counties with county outlines, and a side-by-side combined PNG.
* Included by: gen/fig_walk_il.do, fig_walk_tristate.do, fig_walk_omaha.do, fig_walk_chicago.do.
* Expects: $temp, $figures, $run_tag; ensure_sp_dta; focal() as a list of 5-digit county FIPS.
* Reads:  $temp/tracts_sp.dta; $temp/counties_sp.dta; $temp/final_assignments_<run_tag>_tract.dta
* Writes: $figures/map_<suffix>_hc_<run_tag>.png; $figures/map_<suffix>_leiden_<run_tag>.png; $figures/map_<suffix>_combined_<run_tag>.png
* Notes:  clusters are coloured by a seeded random permutation (seeds 4711 and 4713) over the paired palette, one colour per cluster. Unlike conus_map there is no adjacency-based recolouring, so two bordering clusters can share a hue in a large region. subtitle() and note2() are accepted but not used.

capture program drop region_maps
program define region_maps
    syntax , focal(string) subtitle(string) suffix(string) [note2(string) regionlabel(string)]

    if "`regionlabel'" == "" local regionlabel "region"

    local sp_dta    "$temp/tracts_sp.dta"
    local sp_co_dta "$temp/counties_sp.dta"
    local assign    "$temp/final_assignments_${run_tag}_tract.dta"

    ensure_sp_dta, level(tract)
    if !r(ok) {
        di as text "region_maps (`suffix'): tract shapefile unavailable; skipping."
        exit 0
    }
    ensure_sp_dta, level(county)
    if !r(ok) {
        di as text "region_maps (`suffix'): county shapefile unavailable; skipping."
        exit 0
    }

    capture confirm file "`assign'"
    if _rc {
        di as text "region_maps (`suffix'): `assign' not found (tract level not run); skipping."
        exit 0
    }

    capture frame drop _rm
    capture frame drop _rm_shp
    geoframe create _rm using "`sp_dta'", replace

    frame _rm {
        capture rename GEOID geoid

        gen str5 _ctyfips = substr(geoid, 1, 5)
        gen byte _keep = 0
        foreach c of local focal {
            quietly replace _keep = 1 if _ctyfips == "`c'"
        }
        quietly count if _keep == 1
        display as text "Focal tracts (`regionlabel'): " r(N)
        keep if _keep == 1
        drop _keep _ctyfips

        merge 1:1 geoid using "`assign'", ///
            keep(master match) keepusing(hc_cluster leiden_cluster) generate(_mrg)
        quietly count if _mrg == 1
        if r(N) > 0 {
            display as text "  note: " r(N) " focal tract(s) have no cluster assignment (left uncolored)."
        }
        drop _mrg

        foreach mth in hc leiden {
            local seed = cond("`mth'" == "hc", 4711, 4713)
            egen long _grp = group(`mth'_cluster)
            preserve
            keep _grp
            duplicates drop
            drop if mi(_grp)
            set seed `seed'
            gen double _r = runiform()
            sort _r
            gen long `mth'_color = _n
            keep _grp `mth'_color
            tempfile _perm
            save `_perm'
            restore
            merge m:1 _grp using `_perm', keep(master match) nogenerate
            drop _grp
        }

        quietly summarize hc_color, meanonly
        local n_hc = r(max)
        quietly summarize leiden_color, meanonly
        local n_le = r(max)
        display as text "  Distinct HC clusters in region:      " %5.0f `n_hc'
        display as text "  Distinct Leiden clusters in region:  " %5.0f `n_le'
    }

    capture frame drop _rm_co
    capture frame drop _rm_co_shp
    geoframe create _rm_co using "`sp_co_dta'", replace
    frame _rm_co {
        capture rename GEOID geoid
        gen byte _keep = 0
        foreach c of local focal {
            quietly replace _keep = 1 if geoid == "`c'"
        }
        keep if _keep == 1
        drop _keep
    }

    geoplot (area _rm hc_color, color(paired, n(`n_hc')) lcolor(white) lwidth(vthin)) ///
        (area _rm_co, fcolor(none) lcolor(black) lwidth(medthick)) ///
        , legend(off) name(_hcmap, replace)
    graph export "$figures/map_`suffix'_hc_${run_tag}.png", replace width(2000)
    display as result "  -> map_`suffix'_hc_${run_tag}.png"

    geoplot (area _rm leiden_color, color(paired, n(`n_le')) lcolor(white) lwidth(vthin)) ///
        (area _rm_co, fcolor(none) lcolor(black) lwidth(medthick)) ///
        , legend(off) name(_lemap, replace)
    graph export "$figures/map_`suffix'_leiden_${run_tag}.png", replace width(2000)
    display as result "  -> map_`suffix'_leiden_${run_tag}.png"

    graph combine _hcmap _lemap, ///
        rows(1) ///
        name(_combomap, replace)
    graph export "$figures/map_`suffix'_combined_${run_tag}.png", replace width(4000)
    display as result "  -> map_`suffix'_combined_${run_tag}.png"
    capture graph drop _hcmap _lemap _combomap

    capture frame drop _rm
    capture frame drop _rm_shp
    capture frame drop _rm_co
    capture frame drop _rm_co_shp
end
