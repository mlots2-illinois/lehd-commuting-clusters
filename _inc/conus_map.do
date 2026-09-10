* conus_map.do -- defines conus_map, which draws the HC and Leiden cluster maps of the contiguous United States (or a subset of states) as PDFs with cluster borders.
* Included by: gen/fig_map_national.do, gen/fig_map_regions.do.
* Expects: $temp, $figures, $run_tag, $noncontig_fips, $py, $program; ensure_sp_dta from _inc/ensure_sp_dta.do.
* Reads:  $temp/counties_sp.dta or $temp/tracts_sp.dta; $temp/final_assignments_<run_tag>_<lvl>.dta
* Writes: $figures/map_<suffix>_hc_<run_tag>.pdf; $figures/map_<suffix>_leiden_<run_tag>.pdf; $temp/_gc_nodes_<meth>.csv, _gc_edges_<meth>.csv, _gc_color_<meth>.csv (temporary, erased)
* Notes:  colours are first assigned by a seeded random permutation cycled over 20 tableau colours. A greedy colouring of the cluster adjacency (_py/greedy_color.py) then overwrites them so no two bordering clusters share a colour; the random pass remains as the fallback when the sidecar fails.

capture program drop conus_map
program define conus_map
    syntax , suffix(string) [keepstates(string) level(string)]

    if "`level'" == "" local level "county"
    if !inlist("`level'", "county", "tract") {
        di as error "conus_map: level() must be county or tract."
        exit 198
    }
    if "`level'" == "tract" {
        local sp_dta "$temp/tracts_sp.dta"
        local assign "$temp/final_assignments_${run_tag}_tract.dta"
    }
    else {
        local sp_dta "$temp/counties_sp.dta"
        local assign "$temp/final_assignments_${run_tag}_county.dta"
    }

    ensure_sp_dta, level(`level')
    if !r(ok) {
        di as text "conus_map (`suffix'): `level' shapefile unavailable; skipping."
        exit 0
    }
    capture confirm file "`assign'"
    if _rc {
        di as text "conus_map (`suffix'): `assign' not found (run the `level' cluster pipeline); skipping."
        exit 0
    }

    capture frame drop _nat
    capture frame drop _nat_shp
    geoframe create _nat using "`sp_dta'", replace

    frame _nat {
        capture rename GEOID geoid

        gen str2 _st = substr(geoid, 1, 2)
        gen byte _drop = 0
        foreach s of global noncontig_fips {
            quietly replace _drop = 1 if _st == "`s'"
        }
        if "`keepstates'" != "" {
            gen byte _inkeep = 0
            foreach s of local keepstates {
                quietly replace _inkeep = 1 if _st == "`s'"
            }
            quietly replace _drop = 1 if _inkeep == 0
            drop _inkeep
        }
        quietly drop if _drop == 1
        drop _drop _st

        merge 1:1 geoid using "`assign'", ///
            keep(master match) keepusing(hc_cluster leiden_cluster) nogenerate

        egen long _hc_grp = group(hc_cluster)
        preserve
        keep _hc_grp
        duplicates drop
        drop if mi(_hc_grp)
        set seed 4711
        gen double _r = runiform()
        sort _r
        gen long hc_color = mod(_n - 1, 20) + 1
        keep _hc_grp hc_color
        tempfile _hc_perm
        save `_hc_perm'
        restore
        merge m:1 _hc_grp using `_hc_perm', keep(master match) nogenerate
        drop _hc_grp

        egen long _le_grp = group(leiden_cluster)
        preserve
        keep _le_grp
        duplicates drop
        drop if mi(_le_grp)
        set seed 4713
        gen double _r = runiform()
        sort _r
        gen long leiden_color = mod(_n - 1, 20) + 1
        keep _le_grp leiden_color
        tempfile _le_perm
        save `_le_perm'
        restore
        merge m:1 _le_grp using `_le_perm', keep(master match) nogenerate
        drop _le_grp

        quietly summarize hc_color, meanonly
        local n_hc = r(max)
        quietly summarize leiden_color, meanonly
        local n_le = r(max)
        quietly levelsof hc_cluster, local(_hc_levels)
        local n_hc_clust : word count `_hc_levels'
        quietly levelsof leiden_cluster, local(_le_levels)
        local n_le_clust : word count `_le_levels'
        display as text "  [`suffix'] Distinct HC clusters:     " %6.0f `n_hc_clust' "  (cycled over `n_hc' colors)"
        display as text "  [`suffix'] Distinct Leiden clusters: " %6.0f `n_le_clust' "  (cycled over `n_le' colors)"
    }

    frame _nat: geoframe relink

    tempfile _lk
    frame _nat {
        preserve
        keep _ID hc_cluster leiden_cluster
        save `_lk'
        restore
    }

    foreach meth in hc leiden {
        capture frame drop _`meth'_bord
        capture frame drop _`meth'_bord_shp
        frame _nat: geoframe bshare _`meth'_bord _`meth'_bord_shp, unique nodots replace

        frame _`meth'_bord_shp {
            gen long _ord = _n
            rename _ID  _seg
            rename _ID1 _ID
            merge m:1 _ID using `_lk', keep(master match) ///
                keepusing(`meth'_cluster) nogenerate
            rename `meth'_cluster _cl1
            rename _ID _ID1
            rename _ID2 _ID
            merge m:1 _ID using `_lk', keep(master match) ///
                keepusing(`meth'_cluster) nogenerate
            rename `meth'_cluster _cl2
            rename _ID _ID2
            rename _seg _ID
            sort _ord
            drop _ord
            gen byte _between = (_cl1 != _cl2) & !mi(_cl1) & !mi(_cl2)
        }
        local have_`meth'_ol 1
    }

    foreach meth in hc leiden {
        local _nodes "$temp/_gc_nodes_`meth'.csv"
        local _edges "$temp/_gc_edges_`meth'.csv"
        local _out   "$temp/_gc_color_`meth'.csv"
        capture erase "`_out'"

        capture frame drop _gctmp
        frame copy _nat _gctmp
        frame _gctmp {
            keep `meth'_cluster
            rename `meth'_cluster id
            quietly drop if mi(id)
            quietly duplicates drop
            quietly export delimited id using "`_nodes'", replace
        }
        frame drop _gctmp

        capture frame drop _gctmp
        frame copy _`meth'_bord_shp _gctmp
        frame _gctmp {
            quietly keep if _between == 1
            keep _cl1 _cl2
            rename (_cl1 _cl2) (a b)
            quietly duplicates drop
            quietly export delimited a b using "`_edges'", replace
        }
        frame drop _gctmp

        shell "$py" "$program/_py/greedy_color.py" "`_nodes'" "`_edges'" "`_out'"

        capture confirm file "`_out'"
        if !_rc {
            capture frame drop _gc
            frame create _gc
            frame _gc {
                quietly import delimited using "`_out'", varnames(1) clear
                rename id    `meth'_cluster
                rename color _gccolor
                tempfile _gcdta
                quietly save `_gcdta'
            }
            frame drop _gc
            frame _nat {
                capture drop _gccolor
                merge m:1 `meth'_cluster using "`_gcdta'", ///
                    keep(master match) keepusing(_gccolor) nogenerate
                quietly replace `meth'_color = _gccolor if !mi(_gccolor)
                drop _gccolor
            }
        }
        capture erase "`_nodes'"
        capture erase "`_edges'"
        capture erase "`_out'"
    }
    frame _nat: quietly summarize hc_color, meanonly
    local n_hc = r(max)
    frame _nat: quietly summarize leiden_color, meanonly
    local n_le = r(max)

    local hc_ol_layer ""
    if `have_hc_ol' local hc_ol_layer "(line _hc_bord_shp, ifshp(_between) lcolor(black) lwidth(thin))"
    geoplot (area _nat hc_color, color(tableau, n(`n_hc')) lcolor(white) lwidth(vvthin)) ///
        `hc_ol_layer' ///
        , legend(off) name(_natHC, replace)
    graph export "$figures/map_`suffix'_hc_${run_tag}.pdf", replace
    display as result "  -> map_`suffix'_hc_${run_tag}.pdf"

    local le_ol_layer ""
    if `have_leiden_ol' local le_ol_layer "(line _leiden_bord_shp, ifshp(_between) lcolor(black) lwidth(thin))"
    geoplot (area _nat leiden_color, color(tableau, n(`n_le')) lcolor(white) lwidth(vvthin)) ///
        `le_ol_layer' ///
        , legend(off) name(_natLE, replace)
    graph export "$figures/map_`suffix'_leiden_${run_tag}.pdf", replace
    display as result "  -> map_`suffix'_leiden_${run_tag}.pdf"

    capture graph drop _natHC _natLE
    capture frame drop _nat
    capture frame drop _nat_shp
    capture frame drop _hc_bord
    capture frame drop _hc_bord_shp
    capture frame drop _leiden_bord
    capture frame drop _leiden_bord_shp
end
