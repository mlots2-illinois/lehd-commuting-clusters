clear all
cap log close
local _logdate = subinstr("$S_DATE", " ", "", .)
log using "$log/cluster_04_fuse_${run_tag}_`_logdate'.log", replace text

* cluster_04_fuse.do -- links division-local cluster IDs across neighbouring divisions through their shared overlap units and writes a local-to-global cluster map for HC and for Leiden.
* Called by: _master.do (cluster phase, $run_tag = baseline) and _inc/run_cluster.do (robust phase; $run_tag, $div_tag, $period_pflow_suffix and $hc_input set by the calling robust_* script).
* Reads:  $temp/compiled_assignments_<run_tag>_<lvl>.dta
* Writes: $temp/fusion_map_hc_<run_tag>_<lvl>.dta; $temp/fusion_map_leiden_<run_tag>_<lvl>.dta; $temp/_p404_edges_filt_<run_tag>.dta and $temp/_p404_nodes_<run_tag>.dta (scratch, erased at the end)
* Notes:  Two local clusters in adjacent divisions fuse when they share at least min_shared = 2 overlap units and those units make up at least frac_min = 0.5 of the smaller cluster. The count rule stops one stray unit from fusing two clusters; the fraction rule stops a large cluster from absorbing a small one that merely touches it. Fusion is transitive through the Mata union-find, so a chain of pairwise matches becomes one global cluster.

* Fusion needs both: at least min_shared shared overlap units, and those units at least frac_min of the smaller cluster.
local min_shared 2
local frac_min   0.5

capture mata: mata drop _p404_uf_solve()
mata:
    void _p404_uf_solve(real scalar Nnodes_in, real scalar has_edges)
    {
        real colvector parent, roots
        real matrix edges
        real scalar k, a, b, i, x

        // every node starts as its own root; each edge below unions the roots of its two endpoints
        parent = (1::Nnodes_in)

        if (has_edges) {
            edges = st_data(., ("node_A", "node_B"))
            for (k = 1; k <= rows(edges); k++) {
                a = edges[k, 1]
                // find with path halving: each visited node is re-pointed at its grandparent, so later finds are shorter
                while (parent[a] != a) {
                    parent[a] = parent[parent[a]]
                    a = parent[a]
                }
                b = edges[k, 2]
                while (parent[b] != b) {
                    parent[b] = parent[parent[b]]
                    b = parent[b]
                }
                if (a != b) parent[a] = b
            }
        }

        // a final pass resolves every node to its root so Stata can group on it
        roots = J(Nnodes_in, 1, 0)
        for (i = 1; i <= Nnodes_in; i++) {
            x = i
            while (parent[x] != x) x = parent[x]
            roots[i] = x
        }
        st_matrix("Roots", roots)
    }
end

capture program drop _build_fusion
program define _build_fusion
    syntax , src(string) clusvar(name) outfile(string) min_shared(integer) ///
        frac_min(real)

    di as text _n "=== Fusion: `clusvar' ==="

    use "`src'", clear

    quietly {
        bysort geoid: gen byte _multi = (_N > 1)
        count if _multi == 1
        local n_overlap_rows = r(N)
        drop _multi
    }
    di as text "  Overlap rows detected (geoid in >=2 divs): `n_overlap_rows'"

    if `n_overlap_rows' == 0 {
        di as text "  No overlap units; building identity fusion map."

        use "`src'", clear
        keep div_id `clusvar'
        rename `clusvar' local_cluster
        rename div_id    div
        duplicates drop
        sort div local_cluster
        gen long global_cluster = _n
        rename div div_id

        label variable div_id          "Division ID"
        label variable local_cluster   "Cluster ID local to division"
        label variable global_cluster  "Global cluster ID after fusion"

        compress
        save "`outfile'", replace
        di as result "  -> `outfile' (identity, " _N " clusters)"
        exit
    }

    bysort geoid: gen long _n_div = _N
    keep if _n_div >= 2
    keep geoid div_id `clusvar'
    rename (div_id `clusvar') (div_A local_A)

    tempfile _selfcopy
    preserve
    rename (div_A local_A) (div_B local_B)
    save "`_selfcopy'"
    restore

    * Self-join on the overlap units yields every (cluster in div A, cluster in div B) pair that shares a unit; div_A < div_B keeps each pair once.
    joinby geoid using "`_selfcopy'"
    keep if div_A < div_B

    gen byte _one = 1
    gcollapse (count) n_shared = _one, by(div_A local_A div_B local_B)

    preserve
    use "`src'", clear
    gen byte _one2 = 1
    gcollapse (count) _csize = _one2, by(div_id `clusvar')
    tempfile _csizes
    save `_csizes'
    restore

    rename (div_A local_A) (div_id `clusvar')
    merge m:1 div_id `clusvar' using `_csizes', keep(master match) nogenerate
    rename _csize _csizeA
    rename (div_id `clusvar') (div_A local_A)

    rename (div_B local_B) (div_id `clusvar')
    merge m:1 div_id `clusvar' using `_csizes', keep(master match) nogenerate
    rename _csize _csizeB
    rename (div_id `clusvar') (div_B local_B)

    * The share is taken relative to the smaller cluster, so a small cluster lying wholly inside a large one still passes.
    gen double _frac = n_shared / min(_csizeA, _csizeB)

    quietly count if n_shared >= `min_shared'
    local _n_countrule = r(N)
    keep if n_shared >= `min_shared' & _frac >= `frac_min'
    quietly count
    local _n_kept = r(N)
    local _n_blocked = `_n_countrule' - `_n_kept'
    di as text "  Edges passing count rule (n_shared >= `min_shared'): `_n_countrule'"
    di as text "  Edges blocked by fraction guard (frac < `frac_min'): `_n_blocked'"
    di as text "  Edges to union (count and fraction): `_n_kept'"
    drop _csizeA _csizeB _frac

    local edges_filt "$temp/_p404_edges_filt_${run_tag}.dta"
    local nodes      "$temp/_p404_nodes_${run_tag}.dta"

    save "`edges_filt'", replace

    use "`src'", clear
    keep div_id `clusvar'
    rename `clusvar' local_cluster
    rename div_id    div
    duplicates drop
    sort div local_cluster
    gen long node_id = _n
    local Nnodes = _N
    di as text "  Total cluster nodes: `Nnodes'"
    save "`nodes'", replace

    use "`edges_filt'", clear

    rename (div_A local_A) (div local_cluster)
    merge m:1 div local_cluster using "`nodes'", keep(match) nogenerate
    rename node_id node_A
    rename (div local_cluster) (div_A local_A)

    rename (div_B local_B) (div local_cluster)
    merge m:1 div local_cluster using "`nodes'", keep(match) nogenerate
    rename node_id node_B
    rename (div local_cluster) (div_B local_B)

    keep node_A node_B
    sort node_A node_B

    local Nedges  = _N
    local has_edg = (`Nedges' > 0)
    mata: _p404_uf_solve(`Nnodes', `has_edg')

    di as text "  Reloading nodes file: `nodes'"
    use "`nodes'", clear
    svmat long Roots, names(_root)
    rename _root1 root

    * Root IDs are arbitrary node numbers; group() renumbers them 1..G.
    egen long global_cluster = group(root)
    keep div local_cluster global_cluster
    rename div div_id

    label variable div_id          "Division ID"
    label variable local_cluster   "Cluster ID local to division"
    label variable global_cluster  "Global cluster ID after fusion"

    sort div_id local_cluster
    compress
    save "`outfile'", replace

    quietly summarize global_cluster
    di as result "  -> `outfile' (`Nnodes' local clusters -> " r(max) " global)"
end

foreach lvl of global levels {

    display as text _n "{hline 70}"
    display as text "===== cluster_04 LEVEL: `lvl' ====="
    display as text "{hline 70}"

    local src "$temp/compiled_assignments_${run_tag}_`lvl'.dta"

    _build_fusion,                                                ///
        src("`src'")                                              ///
        clusvar(hc_cluster_local)                                 ///
        outfile("$temp/fusion_map_hc_${run_tag}_`lvl'.dta")       ///
        min_shared(`min_shared') frac_min(`frac_min')

    _build_fusion,                                                ///
        src("`src'")                                              ///
        clusvar(leiden_cluster_local)                             ///
        outfile("$temp/fusion_map_leiden_${run_tag}_`lvl'.dta")   ///
        min_shared(`min_shared') frac_min(`frac_min')
}

capture erase "$temp/_p404_edges_filt_${run_tag}.dta"
capture erase "$temp/_p404_nodes_${run_tag}.dta"

display as result _n "cluster_04 done for $levels."

cap log close
