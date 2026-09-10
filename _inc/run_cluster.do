* run_cluster.do -- runs the eight cluster_* scripts in order for the current namespace ($run_tag and the rest), timing each with _inc/_timeit.do.
* Called by: _master.do (jt01 job-type variant, through `TIME'); included by robust/robust_01, 02, 05, 06, 07 and 08 once per variant.
* Expects: $program, $run_tag and the other namespace globals; $timing_file (from _inc/timing_init.do).
* Notes:  every step passes "nototal"; _master.do times the wrapping robust_* script as a whole, so the TOTAL line counts each variant once.

do "$program/_inc/_timeit.do" "  cluster_01_hc[$run_tag]"          "$program/cluster/cluster_01_hc.do"         "nototal"
do "$program/_inc/_timeit.do" "  cluster_02_leiden[$run_tag]"      "$program/cluster/cluster_02_leiden.do"     "nototal"
do "$program/_inc/_timeit.do" "  cluster_03_compile[$run_tag]"     "$program/cluster/cluster_03_compile.do"    "nototal"
do "$program/_inc/_timeit.do" "  cluster_04_fuse[$run_tag]"        "$program/cluster/cluster_04_fuse.do"       "nototal"
do "$program/_inc/_timeit.do" "  cluster_05_reconcile[$run_tag]"   "$program/cluster/cluster_05_reconcile.do"  "nototal"
do "$program/_inc/_timeit.do" "  cluster_06_harmonize[$run_tag]"   "$program/cluster/cluster_06_harmonize.do"  "nototal"
do "$program/_inc/_timeit.do" "  cluster_07_contiguity[$run_tag]"  "$program/cluster/cluster_07_contiguity.do" "nototal"
do "$program/_inc/_timeit.do" "  cluster_08_absorb_singletons[$run_tag]" "$program/cluster/cluster_08_absorb_singletons.do" "nototal"
