* compare_to_baseline.do -- defines compare_to_baseline, which scores an alternative run's final assignments against the baseline run with NMI and ARI, for HC and Leiden separately.
* Included by: gen/tab_division_sweep.do, gen/tab_jobtype_sweep.do, gen/tab_threshold_sweep.do.
* Expects: $temp; the alt() file passed by the caller; compare_partitions (included here from _inc/compare_partitions.do).
* Reads:  $temp/final_assignments_baseline_<lvl>.dta (cached in frame ctb_base_<lvl> across calls); alt() = $temp/final_assignments_<tag>_<lvl>.dta
* Notes:  units are matched 1:1 on geoid and only matched units are scored, so a variant that drops units is compared on the common set. Only the Stata-side r(NMI) and r(ARI) are used; the AMI sidecar still runs on every call.

include "$program/_inc/compare_partitions.do"

capture program drop compare_to_baseline
program define compare_to_baseline, rclass
    syntax , lvl(string) alt(string)

    return scalar NMI_hc = .
    return scalar ARI_hc = .
    return scalar NMI_le = .
    return scalar ARI_le = .

    local bframe "ctb_base_`lvl'"
    capture confirm frame `bframe'
    if _rc {
        local bfile "$temp/final_assignments_baseline_`lvl'.dta"
        capture confirm file "`bfile'"
        if _rc {
            display as text "compare_to_baseline: no baseline at `lvl' (`bfile')."
            exit
        }
        frame create `bframe'
        frame `bframe' {
            use geoid hc_cluster leiden_cluster using "`bfile'", clear
            rename hc_cluster     hc_b
            rename leiden_cluster le_b
        }
    }

    tempname work
    frame copy `bframe' `work'
    frame `work' {
        merge 1:1 geoid using "`alt'", ///
            keepusing(hc_cluster leiden_cluster) keep(match) nogenerate
        rename hc_cluster     hc_p
        rename leiden_cluster le_p

        compare_partitions, a(hc_b) b(hc_p)
        local nmi_hc = r(NMI)
        local ari_hc = r(ARI)
        compare_partitions, a(le_b) b(le_p)
        local nmi_le = r(NMI)
        local ari_le = r(ARI)
    }
    frame drop `work'

    return scalar NMI_hc = `nmi_hc'
    return scalar ARI_hc = `ari_hc'
    return scalar NMI_le = `nmi_le'
    return scalar ARI_le = `ari_le'
end
