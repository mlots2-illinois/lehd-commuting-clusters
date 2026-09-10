* robust_prologue.do -- stashes the ten baseline namespace globals in $robstash_* so a robust_* script can override them and _inc/robust_epilogue.do can put them back.
* Included by: _master.do (before the jt01 run) and robust/robust_01 ... robust_08 (at their start).
* Expects: the globals set by _inc/baseline_namespaces.do and $levels from _master.do.
* Notes:  the stash lives in globals and the epilogue drops them, so a leftover $robstash_* after a run marks a script that ended before its epilogue.

* Everything a robust_* script may override: run_tag and pflow suffix (all), min_flow/flow_tag (05), div_tag/levels (07), job_type/pool_tag (jt01 in _master.do).
global robstash_run_tag             "$run_tag"
global robstash_period_pflow_suffix "$period_pflow_suffix"
global robstash_job_type            "$job_type"
global robstash_pool_tag            "$pool_tag"
global robstash_min_flow            "$min_flow"
global robstash_flow_tag            "$flow_tag"
global robstash_div_tag             "$div_tag"
global robstash_levels              "$levels"
* No robust_* script changes $raw_lodes or $lodes_version; stashing them keeps prologue and epilogue symmetric over the whole baseline namespace.
global robstash_raw_lodes           "$raw_lodes"
global robstash_lodes_version       "$lodes_version"
