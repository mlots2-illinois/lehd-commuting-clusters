clear all
cap log close
local _logdate = subinstr("$S_DATE", " ", "", .)
log using "$log/cluster_08_absorb_singletons_${run_tag}_`_logdate'.log", replace text

* cluster_08_absorb_singletons.do -- merges every one-unit cluster into the adjacent cluster it commutes with most and overwrites final_assignments in place.
* Called by: _master.do (cluster phase, $run_tag = baseline) and _inc/run_cluster.do (robust phase; $run_tag, $div_tag, $period_pflow_suffix and $hc_input set by the calling robust_* script).
* Reads:  $temp/final_assignments_<run_tag>_<lvl>.dta; $temp/pflow_<lvl>_<period_pflow_suffix>.dta; $temp/<lvl>_neighbors.dta
* Writes: $temp/final_assignments_<run_tag>_<lvl>.dta (overwritten)
* Notes:  The move rule and the maxmoves = 5 freeze mirror cluster_07, which runs first because its moves can strand a unit alone. Absorption only enlarges clusters, so a second pass finds only the islands and frozen units left by the first. The hc_was_fragment flags from cluster_07 are kept beside the new hc_was_singleton flags.

capture program drop _absorb_singletons
program define _absorb_singletons, rclass
    syntax , clusvar(name) flagvar(name) nbfile(string) pflow(string) ///
        movevar(name) maxmoves(integer)

    di as text _n "Absorbing singletons of `clusvar' ..."

    sort geoid
    capture drop _unit_id
    gen long _unit_id = _n

    * A singleton is a cluster of size one; the fragment moves in cluster_07 can leave a cluster with one member.
    bysort `clusvar': gen long _csize = _N
    gen byte _is_single = (_csize == 1)
    quietly count if _is_single == 1
    local n_single = r(N)
    quietly replace `flagvar' = 1 if _is_single == 1

    quietly count if _is_single == 1 & `movevar' < `maxmoves'
    local n_movable = r(N)
    display as text "  Singleton units: `n_single' (movable this pass: `n_movable')"

    local n_absorbed = 0
    if `n_movable' > 0 {
        preserve
        keep geoid _unit_id `clusvar'
        rename geoid     geoid_i
        rename _unit_id  tid_i
        rename `clusvar' cluster_i
        tempfile lt_i
        save `lt_i'
        rename geoid_i   geoid_j
        rename tid_i     tid_j
        rename cluster_i cluster_j
        tempfile lt_j
        save `lt_j'

        * Border edges are relabelled with both endpoints' current clusters on every pass, since clusters change between passes.
        use "`nbfile'", clear
        merge m:1 geoid_i using `lt_i', keep(match) nogenerate
        merge m:1 geoid_j using `lt_j', keep(match) nogenerate
        tempfile edges_all
        save `edges_all'
        restore

        preserve
        keep if _is_single == 1 & `movevar' < `maxmoves'
        keep _unit_id `clusvar'
        rename `clusvar' _own_cluster
        tempfile singles
        save `singles'
        restore

        preserve
        use `edges_all', clear
        keep tid_i tid_j cluster_i cluster_j geoid_i geoid_j

        tempfile B
        save `B'
        keep tid_i cluster_j geoid_i geoid_j
        rename (tid_i cluster_j geoid_i geoid_j) ///
            (_unit_id _nbr_cluster _frag_geoid _nbr_geoid)
        tempfile A
        save `A'
        use `B', clear
        keep tid_j cluster_i geoid_j geoid_i
        rename (tid_j cluster_i geoid_j geoid_i) ///
            (_unit_id _nbr_cluster _frag_geoid _nbr_geoid)
        append using `A'

        joinby _unit_id using `singles', unmatched(none)
        keep if _nbr_cluster != _own_cluster

        tempfile reassign
        quietly count
        if r(N) == 0 {
            di as text "  (all movable singletons are isolated islands; " ///
                "leaving them in place)"
            gen long _new_cluster = .
            keep _unit_id _new_cluster
            save `reassign'
        }
        else {
            * pflow keys pairs as geo_i < geo_j; the singleton and neighbour GEOIDs are ordered the same way before the merge.
            gen str20 _cg_i = cond(_frag_geoid < _nbr_geoid, _frag_geoid, _nbr_geoid)
            gen str20 _cg_j = cond(_frag_geoid < _nbr_geoid, _nbr_geoid, _frag_geoid)
            merge m:1 _cg_i _cg_j using "`pflow'", ///
                keep(master match) keepusing(P_ij) nogenerate
            replace P_ij = 0 if mi(P_ij)

            gen byte _one_e = 1
            gcollapse (sum) flow = P_ij (count) n_edges = _one_e, ///
                by(_unit_id _nbr_cluster)

            * Destination: most commuting, then most shared border, then the lower cluster id.
            gsort _unit_id -flow -n_edges _nbr_cluster
            by _unit_id: keep if _n == 1
            keep _unit_id _nbr_cluster
            rename _nbr_cluster _new_cluster
            save `reassign'
        }
        restore

        merge 1:1 _unit_id using `reassign', keep(master match) nogenerate
        quietly count if !mi(_new_cluster)
        local n_absorbed = r(N)
        quietly replace `clusvar' = _new_cluster if !mi(_new_cluster)
        quietly replace `movevar' = `movevar' + 1 if !mi(_new_cluster)
        drop _new_cluster
        display as text "  Singletons absorbed this pass: `n_absorbed'"
    }

    drop _unit_id _csize _is_single
    return scalar n_singletons = `n_single'
    return scalar n_absorbed   = `n_absorbed'
    return scalar n_frozen     = `n_single' - `n_movable'
end

foreach lvl of global levels {

    display as text _n "{hline 70}"
    display as text "===== cluster_08 LEVEL: `lvl' ====="
    display as text "{hline 70}"

    local nb_file "$temp/`lvl'_neighbors.dta"
    capture confirm file "`nb_file'"
    if _rc {
        display as error "cluster_08: neighbour cache `nb_file' missing; run cluster_07 first. Skipping `lvl'."
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
        display as text "  (no pflow file `pflow_src'; absorption uses border edges only)"
    }

    use "$temp/final_assignments_${run_tag}_`lvl'.dta", clear

    gen byte hc_was_singleton     = 0
    gen byte leiden_was_singleton = 0
    gen int  hc_snmoves     = 0
    gen int  leiden_snmoves = 0

    * Same freeze rule as cluster_07: a unit that has moved maxmoves times stays where it is.
    local maxmoves 5
    foreach cv in hc leiden {
        local iter 0
        local moved 1
        while `moved' > 0 & `iter' < 50 {
            local ++iter
            _absorb_singletons, clusvar(`cv'_cluster) flagvar(`cv'_was_singleton) ///
                nbfile("`nb_file'") pflow("`flow_lookup'") ///
                movevar(`cv'_snmoves) maxmoves(`maxmoves')
            local moved  = r(n_absorbed)
            local frozen = r(n_frozen)
            display as text "  [`lvl' `cv'] pass `iter': absorbed `moved', frozen `frozen'"
        }
        if `moved' > 0 {
            display as error "  [`lvl' `cv'] hit pass cap with `moved' still moving (unexpected with freezing)."
        }
        else {
            display as text "  [`lvl' `cv'] singleton absorption converged in `iter' passes."
        }
    }
    drop hc_snmoves leiden_snmoves

    label variable hc_cluster          "HC cluster (contiguity-enforced, singletons absorbed)"
    label variable leiden_cluster      "Leiden cluster (contiguity-enforced, singletons absorbed)"
    label variable hc_was_singleton     "Unit was ever a size-1 HC cluster"
    label variable leiden_was_singleton "Unit was ever a size-1 Leiden cluster"

    order geoid div_id is_core hc_cluster leiden_cluster ///
        hc_was_fragment leiden_was_fragment hc_was_singleton leiden_was_singleton
    sort geoid

    egen byte _tag_hc = tag(hc_cluster)
    quietly count if _tag_hc
    local hc_n = r(N)
    egen byte _tag_le = tag(leiden_cluster)
    quietly count if _tag_le
    local le_n = r(N)
    drop _tag_hc _tag_le
    display as text _n "Post-absorption cluster counts (`lvl'):"
    display as text "  HC clusters:     " %6.0f `hc_n'
    display as text "  Leiden clusters: " %6.0f `le_n'
    display as text "  Units:           " %6.0f _N

    compress
    * Overwrites cluster_07's file so that every downstream script reads one final_assignments per run_tag.
    save "$temp/final_assignments_${run_tag}_`lvl'.dta", replace
    display as result "cluster_08 (`lvl') done. Overwrote final_assignments_${run_tag}_`lvl'.dta."
}

display as result _n "cluster_08 done for $levels."

cap log close
