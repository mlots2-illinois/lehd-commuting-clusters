* timing_done.do -- appends the TOTAL line and the completion stamp to the run-times file and prints where the logs are.
* Included by: _master.do (last line).
* Expects: $timing_file, $timing_total, $log
* Writes: $timing_file (append)

capture file close _timh
file open _timh using "$timing_file", write append
file write _timh "--------------------------------------------------------------" _n
file write _timh %-46s ("TOTAL") "  " %10.1f ($timing_total) _n
file write _timh "completed $S_DATE $S_TIME" _n
file close _timh

display as result _n "Pipeline complete."
display as result "Logs:      $log"
display as result "Run times: $timing_file (total " %1.0f ($timing_total) " s)"
