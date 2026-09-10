* timing_init.do -- opens the per-run timing file and defines the `TIME' local that _master.do uses to run and time each phase script.
* Included by: _master.do (once, after baseline_namespaces.do).
* Expects: $log, $run_tag, $levels, $program
* Writes: $log/run_times_<run_tag>_<date>.txt (header only; _inc/_timeit.do appends the rows)
* Notes:  `TIME' expands to do "$program/_inc/_timeit.do", so `TIME' "<label>" "<file>" ["nototal"] hands label, script and the optional flag to _timeit.do as arguments 1-3.

local _logdate = subinstr("$S_DATE", " ", "", .)
global timing_file  "$log/run_times_${run_tag}_`_logdate'.txt"
global timing_total 0

capture file close _timh
file open _timh using "$timing_file", write replace
file write _timh "LEHD clustering pipeline: per-program run times" _n
file write _timh "run_tag=$run_tag   levels=$levels   started $S_DATE $S_TIME" _n
file write _timh "--------------------------------------------------------------" _n
file write _timh %-46s ("program") "  " %10s ("seconds") _n
file write _timh "--------------------------------------------------------------" _n
file close _timh

local TIME do "$program/_inc/_timeit.do"
