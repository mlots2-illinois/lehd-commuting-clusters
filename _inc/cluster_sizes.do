* cluster_sizes.do -- defines cluster_sizes, which replaces the data in memory with one row per cluster (size, method) for HC and Leiden.
* Included by: gen/fig_cluster_size_ecdf.do, gen/tab_cluster_size.do.
* Expects: hc_cluster and leiden_cluster in memory (a final_assignments_<run_tag>_<lvl>.dta file).

capture program drop cluster_sizes
program define cluster_sizes

    preserve
    bysort hc_cluster: gen long size = _N
    bysort hc_cluster: keep if _n == 1
    keep size
    gen str8 method = "HC"
    tempfile hc_sizes
    save `hc_sizes'
    restore

    bysort leiden_cluster: gen long size = _N
    bysort leiden_cluster: keep if _n == 1
    keep size
    gen str8 method = "Leiden"
    append using `hc_sizes'
end
