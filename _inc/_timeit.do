* _timeit.do -- runs one phase script with `do`, times it, and appends the label and the elapsed seconds to the run-times file.
* Called by: _master.do through the `TIME' local from _inc/timing_init.do, as `TIME' "<label>" "<file>" ["nototal"]; _inc/run_cluster.do calls it directly.
* Expects: $timing_file and $timing_total (set by _inc/timing_init.do); arguments 1-3 as label, script path, optional flag.
* Writes: $timing_file (append)
* Notes:  a third argument of "nototal" keeps the step's seconds out of $timing_total. The robust variants pass it for their inner cluster_* steps, since _master.do already times each robust_* script as a whole.

local _label  `"`1'"'
local _script `"`2'"'

* Wall-clock from Stata's date and time at one-second resolution, so the seconds column is coarse for short scripts.
local _t0 = clock(c(current_date) + " " + c(current_time), "DMYhms")
do `"`_script'"'
local _t1 = clock(c(current_date) + " " + c(current_time), "DMYhms")

local _sec = (`_t1' - `_t0') / 1000
* "nototal" keeps this step out of the TOTAL line; _inc/run_cluster.do passes it so the robust variants are counted once, at the wrapper.
if "`3'" != "nototal"  global timing_total = $timing_total + `_sec'

* Append one row per script, so a run that dies mid-way still leaves the times of the steps that finished.
capture file close _timh
file open _timh using "$timing_file", write append
file write _timh %-46s (`"`_label'"') "  " %10.1f (`_sec') _n
file close _timh

display as text "  [time] `_label': " %9.1f `_sec' " s"
