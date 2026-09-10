* baseline_namespaces.do -- sets the globals that name the baseline run: JT00 jobs, no minimum flow, the 5-year 2019-2023 pool, LODES v8.4.
* Included by: _master.do (once, before timing_init.do).
* Expects: $raw, $yr_start, $yr_end
* Notes:  every robust_* script overrides some of these and _inc/robust_epilogue.do restores them, so this file is the reference definition of the baseline namespace.

global run_tag             "baseline"
global div_tag             "baseline"
global job_type            "JT00"
global min_flow            0
global lodes_version       "v8.4"
global raw_lodes           "$raw/lodes"
global pool_tag            ""
global flow_tag            ""
global period_pflow_suffix "5yr_${yr_start}_${yr_end}"
