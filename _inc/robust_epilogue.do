* robust_epilogue.do -- restores the ten baseline namespace globals stashed by _inc/robust_prologue.do and drops the stash.
* Included by: _master.do (after the jt01 run) and robust/robust_01 ... robust_08 (at their end).
* Expects: the $robstash_* globals set by _inc/robust_prologue.do.
* Notes:  each robust_* script overrides $run_tag and $period_pflow_suffix plus whichever of $min_flow, $flow_tag, $div_tag, $levels, $job_type, $pool_tag it varies. Restoring all ten unconditionally means no variant can leak into the assess phase, even when a script fails part-way; _master.do therefore runs the epilogue after its `capture noisily` block.

* Restore all ten, not only the ones this variant changed, so the assess phase always starts from the baseline namespace.
global run_tag             "$robstash_run_tag"
global period_pflow_suffix "$robstash_period_pflow_suffix"
global job_type            "$robstash_job_type"
global pool_tag            "$robstash_pool_tag"
* min_flow is numeric in _master.do (global min_flow 0), so it is restored without quotes.
global min_flow            $robstash_min_flow
global flow_tag            "$robstash_flow_tag"
global div_tag             "$robstash_div_tag"
global levels              "$robstash_levels"
global raw_lodes           "$robstash_raw_lodes"
global lodes_version       "$robstash_lodes_version"

* Dropping the stash makes the next prologue start clean instead of stashing an already-overridden value.
macro drop robstash_run_tag robstash_period_pflow_suffix robstash_job_type ///
    robstash_pool_tag robstash_min_flow robstash_flow_tag robstash_div_tag ///
    robstash_levels robstash_raw_lodes robstash_lodes_version
