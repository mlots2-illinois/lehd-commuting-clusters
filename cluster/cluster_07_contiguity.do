clear all
cap log close
local _logdate = subinstr("$S_DATE", " ", "", .)
log using "$log/cluster_07_contiguity_${run_tag}_`_logdate'.log", replace text

* cluster_07_contiguity.do -- makes every HC and Leiden cluster spatially contiguous by moving fragment units to the neighbouring cluster they commute with most.
* Called by: _master.do (cluster phase, $run_tag = baseline) and _inc/run_cluster.do (robust phase; $run_tag, $div_tag, $period_pflow_suffix and $hc_input set by the calling robust_* script).
* Reads:  $temp/harmonized_assignments_<run_tag>_<lvl>.dta; $temp/pflow_<lvl>_<period_pflow_suffix>.dta; $temp/<lvl>_neighbors.dta, built once from $shapefiles/us_<lvl>_2020/US_<lvl>_2020.shp via _py/build_contiguity.py
* Writes: $temp/<lvl>_neighbors.dta (cache; $temp/_p407_contig_edges.csv is scratch); $temp/final_assignments_<run_tag>_<lvl>.dta
* Notes:  Within each cluster we find the connected components over Queen-contiguity edges and keep the largest as the cluster proper; every other component is a fragment. A fragment unit moves to the adjacent cluster with the highest summed P_ij, then the most shared border edges. A unit may move at most maxmoves = 5 times and passes stop at 50, so a unit that cycles between clusters freezes rather than loops. Islands with no neighbour in another cluster stay put. The neighbour cache is shared by every run_tag; deleting it forces a rebuild.

capture mata: mata drop _p407_components()
mata:
    void _p407_components(string scalar edge_frame, real scalar Ntracts_in)
    {
        string scalar saved_frame
        real matrix in_edges, both
        real colvector adj_start, positions, nbrs, cursor, comp, queue
        real scalar k, v, i, cum, start, next_comp, qhead, qtail
        real scalar lo, hi, p, u

        saved_frame = st_framecurrent()
        st_framecurrent(edge_frame)
        if (st_nobs() > 0) {
            in_edges = st_data(., ("tid_i", "tid_j"))
        }
        else {
            in_edges = J(0, 2, .)
        }
        st_framecurrent(saved_frame)

        if (rows(in_edges) > 0) {
            // undirected graph: store each edge in both directions, then build a CSR adjacency list from the sorted rows
            both = in_edges \ in_edges[., (2, 1)]
            both = sort(both, (1, 2))
        }
        else {
            both = J(0, 2, .)
        }

        adj_start = J(Ntracts_in + 1, 1, 0)
        for (k = 1; k <= rows(both); k++) {
            v = both[k, 1]
            adj_start[v] = adj_start[v] + 1
        }
        cum = 0
        positions = J(Ntracts_in + 1, 1, 0)
        for (i = 1; i <= Ntracts_in + 1; i++) {
            positions[i] = cum
            cum = cum + adj_start[i]
        }
        nbrs = J(rows(both), 1, 0)
        cursor = positions
        for (k = 1; k <= rows(both); k++) {
            v = both[k, 1]
            cursor[v] = cursor[v] + 1
            nbrs[cursor[v]] = both[k, 2]
        }

        comp = J(Ntracts_in, 1, 0)
        next_comp = 0
        queue = J(Ntracts_in, 1, 0)
        // breadth-first search from each unlabelled node; comp[] receives one id per connected component
        for (start = 1; start <= Ntracts_in; start++) {
            if (comp[start] != 0) continue
            next_comp = next_comp + 1
            comp[start] = next_comp
            qhead = 1; qtail = 1
            queue[1] = start
            while (qhead <= qtail) {
                v = queue[qhead]; qhead = qhead + 1
                lo = positions[v] + 1
                hi = positions[v + 1]
                for (p = lo; p <= hi; p++) {
                    u = nbrs[p]
                    if (comp[u] == 0) {
                        comp[u] = next_comp
                        qtail = qtail + 1
                        queue[qtail] = u
                    }
                }
            }
        }
        st_store(., "_comp", comp)
    }
end

capture program drop _enforce_contig
program define _enforce_contig, rclass
    syntax , clusvar(name) flagvar(name) nbfile(string) pflow(string) ///
        movevar(name) maxmoves(integer)

    di as text _n "Enforcing contiguity on `clusvar' ..."

    sort geoid
    gen long _tract_id = _n

    preserve
    keep geoid _tract_id `clusvar'
    rename geoid     geoid_i
    rename _tract_id tid_i
    rename `clusvar' cluster_i
    tempfile lt_i
    save `lt_i'

    rename geoid_i   geoid_j
    rename tid_i     tid_j
    rename cluster_i cluster_j
    tempfile lt_j
    save `lt_j'

    use "`nbfile'", clear
    merge m:1 geoid_i using `lt_i', keep(match) nogenerate
    merge m:1 geoid_j using `lt_j', keep(match) nogenerate

    gen byte same_cluster = (cluster_i == cluster_j)

    tempfile edges_all
    save `edges_all'
    restore

    preserve
    * Components are found over same-cluster edges only, so a cluster cut in two by another cluster shows as two components.
    use `edges_all' if same_cluster == 1, clear
    keep tid_i tid_j
    tempfile in_edges
    save `in_edges'
    restore

    quietly count
    local Ntracts = r(N)

    capture drop _comp
    gen long _comp = .

    capture frame drop _p407ef
    frame create _p407ef
    frame _p407ef: use `in_edges', clear

    mata: _p407_components("_p407ef", `Ntracts')

    frame drop _p407ef

    bysort `clusvar' _comp: gen long _comp_size = _N
    * The largest component is the cluster proper; ties go to the higher component id.
    bysort `clusvar' (_comp_size _comp): ///
        gen long _max_comp = _comp[_N]
    gen byte _is_fragment = (_comp != _max_comp)
    quietly replace `flagvar' = 1 if _is_fragment == 1
    quietly count if _is_fragment == 1
    local n_frag = r(N)
    quietly count if _is_fragment == 1 & `movevar' < `maxmoves'
    local n_movable = r(N)
    display as text "  Fragment units: `n_frag' (movable this pass: `n_movable')"

    local n_reassigned = 0
    if `n_movable' > 0 {
        preserve
        keep if _is_fragment == 1 & `movevar' < `maxmoves'
        keep _tract_id `clusvar'
        rename `clusvar' _own_cluster
        tempfile fragments
        save `fragments'
        restore

        preserve
        use `edges_all', clear
        keep tid_i tid_j cluster_i cluster_j geoid_i geoid_j

        tempfile B
        save `B'

        keep tid_i cluster_j geoid_i geoid_j
        rename (tid_i cluster_j geoid_i geoid_j) ///
            (_tract_id _nbr_cluster _frag_geoid _nbr_geoid)
        tempfile A
        save `A'

        use `B', clear
        keep tid_j cluster_i geoid_j geoid_i
        rename (tid_j cluster_i geoid_j geoid_i) ///
            (_tract_id _nbr_cluster _frag_geoid _nbr_geoid)
        append using `A'

        joinby _tract_id using `fragments', unmatched(none)
        keep if _nbr_cluster != _own_cluster

        tempfile reassign
        quietly count
        if r(N) == 0 {
            di as text "  (all fragment units are isolated islands; " ///
                "leaving them in place)"
            gen long _new_cluster = .
            keep _tract_id _new_cluster
            save `reassign'
        }
        else {
            * pflow keys pairs as geo_i < geo_j; the fragment and neighbour GEOIDs are ordered the same way before the merge.
            gen str20 _cg_i = cond(_frag_geoid < _nbr_geoid, _frag_geoid, _nbr_geoid)
            gen str20 _cg_j = cond(_frag_geoid < _nbr_geoid, _nbr_geoid, _frag_geoid)
            merge m:1 _cg_i _cg_j using "`pflow'", ///
                keep(master match) keepusing(P_ij) nogenerate
            replace P_ij = 0 if mi(P_ij)

            gen byte _one_e = 1
            gcollapse (sum) flow = P_ij (count) n_edges = _one_e, ///
                by(_tract_id _nbr_cluster)

            * Destination: most commuting, then most shared border, then the lower cluster id.
            gsort _tract_id -flow -n_edges _nbr_cluster
            by _tract_id: keep if _n == 1
            keep _tract_id _nbr_cluster
            rename _nbr_cluster _new_cluster
            save `reassign'
        }
        restore

        merge 1:1 _tract_id using `reassign', keep(master match) nogenerate
        quietly count if !mi(_new_cluster)
        local n_reassigned = r(N)
        quietly replace `clusvar' = _new_cluster if !mi(_new_cluster)
        quietly replace `movevar' = `movevar' + 1 if !mi(_new_cluster)
        drop _new_cluster
        display as text "  Fragment units reassigned this pass: `n_reassigned'"
    }

    drop _tract_id _comp _comp_size _max_comp _is_fragment
    return scalar n_fragments  = `n_frag'
    return scalar n_reassigned = `n_reassigned'
    return scalar n_frozen     = `n_frag' - `n_movable'
end

capture program drop _ensure_neighbors
program define _ensure_neighbors
    syntax , shp_file(string) nb_file(string)

    capture confirm file "`nb_file'"
    if !_rc exit

    capture confirm file "`shp_file'"
    if _rc {
        di as error "Shapefile not found: `shp_file'"
        di as error "Cannot enforce contiguity at this level; skipping."
        exit 601
    }

    display as text "Building Queen contiguity neighbors (one-time) -> `nb_file' ..."

    local edges_csv "$temp/_p407_contig_edges.csv"
    display as text "  Calling: $py build_contiguity.py ..."
    shell "$py" "$program/_py/build_contiguity.py" "`shp_file'" "`edges_csv'" GEOID

    capture confirm file "`edges_csv'"
    if _rc {
        di as error "Python contiguity helper did not produce `edges_csv'."
        di as error "Make sure geopandas + libpysal are installed in $py."
        exit 601
    }

    import delimited using "`edges_csv'", varnames(1) ///
        stringcols(1 2) clear
    sort geoid_i geoid_j
    compress
    save "`nb_file'", replace
    erase "`edges_csv'"
    display as text "  Cached " _N " Queen edges -> `nb_file'"
end

foreach lvl of global levels {

    display as text _n "{hline 70}"
    display as text "===== cluster_07 LEVEL: `lvl' ====="
    display as text "{hline 70}"

    if "`lvl'" == "tract" {
        local shp_dir   "$shapefiles/us_tract_2020"
        local shp_file  "`shp_dir'/US_tract_2020.shp"
        local nb_file   "$temp/tract_neighbors.dta"
    }
    else {
        local shp_dir   "$shapefiles/us_county_2020"
        local shp_file  "`shp_dir'/US_county_2020.shp"
        local nb_file   "$temp/county_neighbors.dta"
    }

    capture noisily _ensure_neighbors,              ///
        shp_file("`shp_file'") nb_file("`nb_file'")
    if _rc {
        display as error "Skipping cluster_07 contiguity for `lvl' (rc=" _rc ")."
        continue
    }

    local pflow_src "$temp/pflow_`lvl'_${period_pflow_suffix}.dta"
    tempfile flow_lookup
    capture confirm file "`pflow_src'"
    if !_rc {
        preserve
        use geo_i geo_j P_ij using "`pflow_src'", clear
        rename (geo_i geo_j) (_cg_i _cg_j)
        recast str20 _cg_i _cg_j
        save "`flow_lookup'"
        restore
    }
    else {
        preserve
        clear
        set obs 0
        gen str20 _cg_i = ""
        gen str20 _cg_j = ""
        gen double P_ij = .
        save "`flow_lookup'"
        restore
        display as text "  (no pflow file `pflow_src'; reassignment uses border edges only)"
    }

    use "$temp/harmonized_assignments_${run_tag}_`lvl'.dta", clear
    gen byte hc_was_fragment     = 0
    gen byte leiden_was_fragment = 0
    gen int  hc_nmoves     = 0
    gen int  leiden_nmoves = 0

    * A unit that has moved maxmoves times is frozen; without the cap two clusters can trade a unit back and forth forever.
    local maxmoves 5
    foreach cv in hc leiden {
        local iter 0
        local reas 1
        while `reas' > 0 & `iter' < 50 {
            local ++iter
            _enforce_contig, clusvar(`cv'_cluster) flagvar(`cv'_was_fragment) ///
                nbfile("`nb_file'") pflow("`flow_lookup'") ///
                movevar(`cv'_nmoves) maxmoves(`maxmoves')
            local reas   = r(n_reassigned)
            local frozen = r(n_frozen)
            display as text "  [`lvl' `cv'] pass `iter': reassigned `reas', frozen `frozen'"
        }
        if `reas' > 0 {
            display as error "  [`lvl' `cv'] hit pass cap with `reas' still moving (unexpected with freezing)."
        }
        else {
            display as text "  [`lvl' `cv'] contiguity converged in `iter' passes."
        }
    }
    drop hc_nmoves leiden_nmoves

    order geoid div_id is_core hc_cluster leiden_cluster hc_was_fragment leiden_was_fragment
    sort geoid

    label variable hc_cluster          "HC cluster (contiguity-enforced)"
    label variable leiden_cluster      "Leiden cluster (contiguity-enforced)"
    label variable hc_was_fragment     "Unit was ever part of a non-contiguous HC fragment"
    label variable leiden_was_fragment "Unit was ever part of a non-contiguous Leiden fragment"

    compress
    save "$temp/final_assignments_${run_tag}_`lvl'.dta", replace
    display as result "cluster_07 (`lvl') done. Saved final_assignments_${run_tag}_`lvl'.dta."
}

display as result _n "cluster_07 done for $levels."

cap log close
