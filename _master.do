* _master.do -- entry point for the LODES commuting-cluster pipeline; sets paths and parameters, then runs the phases in order.
* Run as:  stata-mp -b do _master.do  with LEHD_PROJ, LEHD_PROG and LEHD_PYENV exported (see README.md).
* Phases:  data -> build -> calib -> cluster -> robust -> assess -> generate -> export, each behind a do_* switch below.
* Writes:  $output/logs/run_times_<run_tag>_<date>.txt and, through the phase scripts, everything under $temp and $output.
* Notes:   Every intermediate carries $run_tag in its name; the robust phase reruns the cluster phase under other tags and restores the baseline afterwards.

clear all
set more off
version 18

local _missing ""
foreach _cmd in gcollapse gegen gcontract geoframe geoplot spshape2dta {
    capture which `_cmd'
    if _rc local _missing "`_missing' `_cmd'"
}

if "`_missing'" != "" {
    display as error _n "{hline 70}"
    display as error "MISSING STATA DEPENDENCIES:`_missing'"
    exit 111
}

global proj : environment LEHD_PROJ
if "$proj" == "" global proj "/path/to/project"

global program : environment LEHD_PROG
if "$program" == "" global program "/path/to/lehd-commuting-clusters"

global raw       "$proj/03_data/01_raw"
global shapefiles "$proj/03_data/00_shapefiles"
global temp      "$proj/03_data/02_processed"

global output : environment LEHD_OUT
if "$output" == "" global output "$proj/05_output"

global log       "$output/logs"
global figures   "$output/figures"
global tables    "$output/tables"

capture mkdir "$temp"
capture mkdir "$output"
capture mkdir "$log"
capture mkdir "$figures"
capture mkdir "$tables"

global pyenv : environment LEHD_PYENV
if "$pyenv" == "" global pyenv "/path/to/python-env"

global py : environment LEHD_PY
if "$py" == "" global py "$pyenv/bin/python"

global divfmt    "%05.0f"

global yr_start    2019
global yr_end      2023

global seed        202600331

global noncontig_fips "02 15 60 66 69 72 78"

include "$program/_inc/div_specs.do"

numlist "0(0.05)1"
global alpha_sc_grid "`r(numlist)'"

global levels      "tract county"
global job_types   "JT00 JT01"

global do_data     1
global do_build    1
global do_calib    1
global do_cluster  1
global do_robust   1
global do_robust_window     1
global do_robust_jobtype    1
global do_robust_threshold  1
global do_robust_params     1
global do_robust_divisions  1
global do_robust_alpha_full 1

global do_assess   1
global do_generate 1
global do_export   1
global do_sync_manuscript 0

include "$program/_inc/baseline_namespaces.do"

include "$program/_inc/timing_init.do"

if $do_data {
    `TIME' "data_01_state_lookup" "$program/data/data_01_state_lookup.do"

    foreach jt of global job_types {
        global job_type "`jt'"
        global pool_tag = cond("`jt'" == "JT00", "", "_" + lower("`jt'"))
        `TIME' "data_02_create_files[`jt']" "$program/data/data_02_create_files.do"
        `TIME' "data_03_pool_multi[`jt']"   "$program/data/data_03_pool_multi.do"
    }

    global job_type "JT00"
    global pool_tag ""

    `TIME' "data_04_cbsa_crosswalk" "$program/data/data_04_cbsa_crosswalk.do"
}

if $do_build {

    foreach jt of global job_types {
        global job_type "`jt'"
        global pool_tag = cond("`jt'" == "JT00", "", "_" + lower("`jt'"))
        `TIME' "build_01_pflows[`jt']" "$program/build/build_01_pflows.do"
    }
    global job_type "JT00"
    global pool_tag ""

    foreach m of numlist 2 3 5 10 20 {
        global min_flow `m'
        global flow_tag "_mf`m'"
        `TIME' "build_01_pflows[mf`m']" "$program/build/build_01_pflows.do"
    }
    global min_flow 0
    global flow_tag ""

    `TIME' "build_02_centroids" "$program/build/build_02_centroids.do"

    foreach lvl of global levels {
        foreach s of global div_specs_`lvl' {
            gettoken dtag  rest : s
            gettoken bands rest : rest
            gettoken W     S    : rest

            global div_tag    "`dtag'"
            global div_bands  `bands'
            global div_levels "`lvl'"
            if "`lvl'" == "tract" {
                global div_W `W'
                global div_S `S'
            }
            else {
                global div_W_county `W'
                global div_S_county `S'
            }

            display as result _n "{hline 70}"
            display as result ///
                "build_03: `lvl' division geometry '`dtag'' (bands=`bands' W=`W' S=`S')"
            display as result "{hline 70}"
            `TIME' "build_03_divisions[`lvl'/`dtag']" "$program/build/build_03_divisions.do"
        }
    }
    global div_tag    "baseline"
    global div_levels ""

    `TIME' "data_04_cbsa_crosswalk[tract]" "$program/data/data_04_cbsa_crosswalk.do"
}

if $do_calib {

    global n_sample_divs 25

    global do_calib_sweeps 1

    if $do_calib_sweeps {

        numlist "0(0.05)1"
        global alpha_grid       "`r(numlist)'"
        numlist "10(10)300"
        global cap_grid_tract   "`r(numlist)'"
        numlist "50(10)2000"
        global cap_grid_county  "`r(numlist)'"
        global disp_cap_rows_tract   "10 50 100 150 200 300"
        global disp_cap_rows_county  "50 500 1000 1320 1600 2000"
        global disp_alpha_rows       "0 0.20 0.40 0.60 0.80 0.90 1"

        * The scripts run 01, 03, 02. The file numbers follow the order in
        * which the paper presents the sweeps; neither 02 nor 03 reads the
        * other's output, so the run order is free.
        `TIME' "calib_01_alpha_cap" "$program/calib/calib_01_alpha_cap.do"

        numlist "5000(1000)15000 17500(2500)40000", sort
        global k_grid_tract "`r(numlist)'"

        numlist "500(10)1000"
        global k_grid_county "`r(numlist)'"

        `TIME' "calib_03_hc_k" "$program/calib/calib_03_hc_k.do"

        numlist "0.100(0.005)0.200 0.171", sort
        local clean_list = subinstr("`r(numlist)'", " ", ",", .)
        global resolutions_tract  "`clean_list'"

        numlist "0.150(0.002)0.300", sort
        local clean_list = subinstr("`r(numlist)'", " ", ",", .)
        global resolutions_county "`clean_list'"

        global leiden_sweep_step 0.010

        global leiden_nseeds_tract  10
        global leiden_nseeds_county 10

        `TIME' "calib_02_resolution" "$program/calib/calib_02_resolution.do"
    }

    `TIME' "calib_04_pick_params" "$program/calib/calib_04_pick_params.do"

    * cluster-validity evidence for the chosen cut: bootstrap stability
    * (Hennig 2007) and the geography-only / permutation-null benchmarks.
    * Runs after calib_04 because it reads the chosen k.
    global do_calib_validity 1
    global hs_nboot 50
    global nt_nperm 99
    if $do_calib_validity {
        `TIME' "calib_05_validity" "$program/calib/calib_05_validity.do"
    }
}

global hc_input  "$temp/hc_input"

if $do_cluster {
    `TIME' "cluster_01_hc[$run_tag]"                    "$program/cluster/cluster_01_hc.do"
    `TIME' "cluster_02_leiden[$run_tag]"                "$program/cluster/cluster_02_leiden.do"
    `TIME' "cluster_03_compile[$run_tag]"               "$program/cluster/cluster_03_compile.do"
    `TIME' "cluster_04_fuse[$run_tag]"                  "$program/cluster/cluster_04_fuse.do"
    `TIME' "cluster_05_reconcile[$run_tag]"             "$program/cluster/cluster_05_reconcile.do"
    `TIME' "cluster_06_harmonize[$run_tag]"             "$program/cluster/cluster_06_harmonize.do"
    `TIME' "cluster_07_contiguity[$run_tag]"            "$program/cluster/cluster_07_contiguity.do"
    `TIME' "cluster_08_absorb_singletons[$run_tag]"     "$program/cluster/cluster_08_absorb_singletons.do"
}

if $do_robust {
    if $do_robust_window {
        `TIME' "robust_01_window_3yr"               "$program/robust/robust_01_window_3yr.do"
        `TIME' "robust_02_window_1yr"               "$program/robust/robust_02_window_1yr.do"
        `TIME' "robust_03_window_baseline_ref"      "$program/robust/robust_03_window_baseline_ref.do"
    }

    if $do_robust_threshold  `TIME' "robust_05_threshold" "$program/robust/robust_05_threshold.do"
    if $do_robust_params     `TIME' "robust_06_params"    "$program/robust/robust_06_params.do"
    if $do_robust_divisions  `TIME' "robust_07_divisions" "$program/robust/robust_07_divisions.do"
    if $do_robust_alpha_full `TIME' "robust_08_alpha_full" "$program/robust/robust_08_alpha_full.do"

    if $do_robust_jobtype {
        include "$program/_inc/robust_prologue.do"
        global job_type            "JT01"
        global pool_tag            "_jt01"
        global run_tag             "jt01"
        global period_pflow_suffix "5yr_${yr_start}_${yr_end}_jt01"
        capture noisily `TIME' "run_cluster[jt01]" "$program/_inc/run_cluster.do"
        local _jt01_rc = _rc
        include "$program/_inc/robust_epilogue.do"
        if `_jt01_rc' ///
            display as error "robust[jt01]: failed (rc=`_jt01_rc'); namespaces restored."
    }
}

if $do_assess {
    `TIME' "assess_01_elbow"                 "$program/assess/assess_01_elbow.do"
    `TIME' "assess_02_agreement"             "$program/assess/assess_02_agreement.do"
    `TIME' "assess_03_temporal"              "$program/assess/assess_03_temporal.do"
    `TIME' "assess_04_sensitivity"           "$program/assess/assess_04_sensitivity.do"
    `TIME' "assess_05_cbsa"                  "$program/assess/assess_05_cbsa.do"
    `TIME' "assess_06_czones"                "$program/assess/assess_06_czones.do"
    `TIME' "assess_07_partition_quality"     "$program/assess/assess_07_partition_quality.do"
    `TIME' "assess_08_partition_quality_alpha" "$program/assess/assess_08_partition_quality_alpha.do"
    `TIME' "assess_09_maup"                  "$program/assess/assess_09_maup.do"
    `TIME' "assess_10_maup_paramsweep"       "$program/assess/assess_10_maup_paramsweep.do"
    `TIME' "assess_11_robinson"             "$program/assess/assess_11_robinson.do"
    `TIME' "assess_12_random_null"          "$program/assess/assess_12_random_null.do"
}

if $do_generate {

    `TIME' "tab_chosen_params"           "$program/gen/tab_chosen_params.do"
    `TIME' "tab_cluster_size"            "$program/gen/tab_cluster_size.do"
    `TIME' "tab_leiden_sweep"            "$program/gen/tab_leiden_sweep.do"
    `TIME' "tab_temporal_stability"      "$program/gen/tab_temporal_stability.do"
    `TIME' "tab_temporal_yoy"            "$program/gen/tab_temporal_yoy.do"
    `TIME' "tab_alpha_sweep"             "$program/gen/tab_alpha_sweep.do"
    `TIME' "tab_dbar_sweep"              "$program/gen/tab_dbar_sweep.do"
    `TIME' "tab_hc_k_sensitivity"        "$program/gen/tab_hc_k_sensitivity.do"
    `TIME' "tab_leiden_gamma_sensitivity" "$program/gen/tab_leiden_gamma_sensitivity.do"
    `TIME' "tab_jobtype_sweep"           "$program/gen/tab_jobtype_sweep.do"
    `TIME' "tab_threshold_sweep"         "$program/gen/tab_threshold_sweep.do"
    `TIME' "tab_division_sweep"          "$program/gen/tab_division_sweep.do"

    `TIME' "fig_snakescan"               "$program/gen/fig_snakescan.do"
    `TIME' "fig_admin_places"            "$program/gen/fig_admin_places.do"
    `TIME' "fig_calib_dispersion"        "$program/gen/fig_calib_dispersion.do"
    `TIME' "fig_gamma_fine"              "$program/gen/fig_gamma_fine.do"
    `TIME' "fig_elbow_hc_clusters"       "$program/gen/fig_elbow_hc_clusters.do"
    `TIME' "fig_elbow_leiden_clusters"   "$program/gen/fig_elbow_leiden_clusters.do"
    `TIME' "fig_elbow_overlay"           "$program/gen/fig_elbow_overlay.do"
    `TIME' "fig_elbow_clusters"          "$program/gen/fig_elbow_clusters.do"
    `TIME' "fig_leiden_sweep_diagnostic" "$program/gen/fig_leiden_sweep_diagnostic.do"
    `TIME' "fig_cluster_size_ecdf"       "$program/gen/fig_cluster_size_ecdf.do"
    `TIME' "fig_sensitivity_dotplot"     "$program/gen/fig_sensitivity_dotplot.do"
    `TIME' "fig_maup_sensitivity"        "$program/gen/fig_maup_sensitivity.do"
    `TIME' "fig_maup_paramsweep"         "$program/gen/fig_maup_paramsweep.do"
    `TIME' "fig_walk_il"                 "$program/gen/fig_walk_il.do"
    `TIME' "fig_walk_tristate"           "$program/gen/fig_walk_tristate.do"
    `TIME' "fig_walk_omaha"              "$program/gen/fig_walk_omaha.do"
    `TIME' "fig_walk_chicago"            "$program/gen/fig_walk_chicago.do"
    `TIME' "fig_map_national"            "$program/gen/fig_map_national.do"
    `TIME' "fig_map_regions"             "$program/gen/fig_map_regions.do"
    `TIME' "fig_map_states"              "$program/gen/fig_map_states.do"
    `TIME' "fig_flow_map"                "$program/gen/fig_flow_map.do"

    `TIME' "report_01_descriptive_table" "$program/report/report_01_descriptive_table.do"
    `TIME' "report_02_dissim_example"    "$program/report/report_02_dissim_example.do"
    `TIME' "report_04_validation"        "$program/report/report_04_validation.do"
}

if $do_export {
    `TIME' "report_10_export" "$program/export/report_10_export.do"
}

include "$program/_inc/timing_done.do"
