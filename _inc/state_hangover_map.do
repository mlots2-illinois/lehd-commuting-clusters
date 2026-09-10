* state_hangover_map.do -- defines state_hangover_map, which maps every tract-level cluster that touches a given state, so clusters crossing the state line are drawn whole.
* Included by: gen/fig_map_states.do.
* Expects: $temp, $figures, $run_tag, $py, $program; ensure_sp_dta.
* Reads:  $temp/tracts_sp.dta; $temp/states_sp.dta (state outline, optional); $temp/final_assignments_<run_tag>_tract.dta
* Writes: $figures/map_state_<suffix>_<meth>_<run_tag>.pdf for meth in hc, leiden; $temp/_gcs_nodes.csv, _gcs_edges.csv, _gcs_color.csv (temporary, erased)
* Notes:  colouring follows conus_map: a seeded random cycle over 20 colours, then greedy recolouring of the cluster adjacency through _py/greedy_color.py.

capture program drop state_hangover_map
program define state_hangover_map
    syntax , state(string) suffix(string)

    local assign "$temp/final_assignments_${run_tag}_tract.dta"

    ensure_sp_dta, level(tract)
    if !r(ok) {
        di as text "state_hangover_map (`suffix'): tract shapefile unavailable; skipping."
        exit 0
    }
    capture confirm file "`assign'"
    if _rc {
        di as text "state_hangover_map (`suffix'): `assign' not found (run the tract pipeline); skipping."
        exit 0
    }
    ensure_sp_dta, level(state)
    local have_state = r(ok)

    if `have_state' {
        capture frame drop _stbnd
        capture frame drop _stbnd_shp
        geoframe create _stbnd using "$temp/states_sp.dta", replace
        frame _stbnd: keep if GEOID == "`state'"
    }

    foreach meth in hc leiden {

        capture frame drop _nat
        capture frame drop _nat_shp
        geoframe create _nat using "$temp/tracts_sp.dta", replace
        frame _nat {
            capture rename GEOID geoid
            merge 1:1 geoid using "`assign'", ///
                keep(master match) keepusing(`meth'_cluster) nogenerate
            gen byte _instate = (substr(geoid, 1, 2) == "`state'")
            egen byte _touch = max(_instate), by(`meth'_cluster)
            quietly drop if _touch != 1 | mi(`meth'_cluster)
            drop _instate _touch

            egen long _grp = group(`meth'_cluster)
            preserve
            keep _grp
            duplicates drop
            drop if mi(_grp)
            set seed 4719
            gen double _r = runiform()
            sort _r
            gen long `meth'_color = mod(_n - 1, 20) + 1
            keep _grp `meth'_color
            tempfile _perm
            save `_perm'
            restore
            merge m:1 _grp using `_perm', keep(master match) nogenerate
            drop _grp
            quietly levelsof `meth'_cluster, local(_lv)
            display as text "  [`suffix' `meth'] clusters touching state: " ///
                `:word count `_lv'' "  (tracts " _N ")"
        }

        frame _nat: geoframe relink

        tempfile _lk
        frame _nat {
            preserve
            keep _ID `meth'_cluster
            save `_lk'
            restore
        }
        capture frame drop _bord
        capture frame drop _bord_shp
        frame _nat: geoframe bshare _bord _bord_shp, unique nodots replace
        frame _bord_shp {
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

        local _nodes "$temp/_gcs_nodes.csv"
        local _edges "$temp/_gcs_edges.csv"
        local _out   "$temp/_gcs_color.csv"
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
        frame copy _bord_shp _gctmp
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

        frame _nat: quietly summarize `meth'_color, meanonly
        local ncol = r(max)

        local ol ""
        frame _bord_shp: quietly count if _between == 1
        if r(N) > 0 ///
            local ol "(line _bord_shp, ifshp(_between) lcolor(black) lwidth(thin))"
        local stlayer ""
        if `have_state' ///
            local stlayer "(area _stbnd, fcolor(none) lcolor(black) lwidth(medthin))"

        geoplot (area _nat `meth'_color, color(tableau, n(`ncol')) ///
            lcolor(white) lwidth(vvthin)) ///
            `ol' `stlayer' ///
            , legend(off) name(_stmap, replace)
        graph export "$figures/map_state_`suffix'_`meth'_${run_tag}.pdf", replace
        display as result "  -> map_state_`suffix'_`meth'_${run_tag}.pdf"

        capture graph drop _stmap
        capture frame drop _nat
        capture frame drop _nat_shp
        capture frame drop _bord
        capture frame drop _bord_shp
    }
    capture frame drop _stbnd
    capture frame drop _stbnd_shp
end
