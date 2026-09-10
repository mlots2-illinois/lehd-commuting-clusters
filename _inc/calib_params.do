* calib_params.do -- hand-set fallback values of the four calibrated parameters per level: alpha, cap_km, hc_k_national and leiden_res.
* Included by: calib/calib_01_alpha_cap.do, calib_02_resolution.do, calib_03_hc_k.do, calib_04_pick_params.do, calib_05_validity.do, assess/assess_08_partition_quality_alpha.do, gen/fig_leiden_sweep_diagnostic.do, _inc/sens_grid.do.
* Expects: nothing; defines the locals `tract_*' and `county_*' in the including scope.
* Notes:  where $temp/chosen_params_<lvl>.dta exists (written by calib_04) the callers read it and these locals are only the fallback. calib_02 prints a reminder to set tract_leiden_res and county_leiden_res here by hand after the sweep.

local tract_alpha            0.8
local tract_cap_km           100

local tract_hc_k_national    9000
local tract_leiden_res       0.171

local county_alpha            0.8
local county_cap_km           1320

local county_hc_k_national    700
local county_leiden_res       0.226
