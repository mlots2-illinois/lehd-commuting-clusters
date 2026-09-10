* div_specs.do -- lists the division geometries build_03_divisions.do builds, per level: tag, number of latitude bands, window W and stride S.
* Included by: _master.do (once, before the build phase).
* Expects: nothing.
* Notes:  the tract list holds the baseline (30 bands, W = 3000, S = 2000) and fourteen variants that robust_07_divisions.do re-clusters. County has the baseline only, with W = S = 5000, which puts all counties in one window.

global div_specs_tract `" "baseline 30 3000 2000" "'
global div_specs_tract `"$div_specs_tract "Wnarrow 30 2000 1000" "Wmid 30 2500 1500" "'
global div_specs_tract `"$div_specs_tract "Wwide 30 4000 3000" "Wxwide 30 5000 4000" "'
global div_specs_tract `"$div_specs_tract "ovxsmall 30 3000 2750" "ovsmall 30 3000 2500" "'
global div_specs_tract `"$div_specs_tract "ovlarge 30 3000 1500" "ovxlarge 30 3000 1000" "'
global div_specs_tract `"$div_specs_tract "bands15 15 3000 2000" "bands20 20 3000 2000" "'
global div_specs_tract `"$div_specs_tract "bands45 45 3000 2000" "bands60 60 3000 2000" "'
global div_specs_tract `"$div_specs_tract "smallboth 30 2000 1500" "bigboth 30 4000 2000" "'
global div_specs_county `" "baseline 30 5000 5000" "'
