# lehd-commuting-clusters

Replication code for delineating commuting-based labor markets in the contiguous United States from LEHD Origin-Destination Employment Statistics (LODES). We pool five years of LODES version 8 origin-destination (OD) flows, 2019 through 2023, for the 48 contiguous states and the District of Columbia. We then partition census tracts and counties with two algorithms run on the same flow-distance dissimilarity: agglomerative hierarchical clustering (HC) and Leiden community detection under the constant Potts model (CPM). The two partitions are fused, reconciled across geographies, made spatially contiguous, and exported as a tract-level and a county-level assignment file. The remainder of the pipeline calibrates the parameters, reruns the delineation under alternative windows, thresholds, parameters, and division geometries, and assesses the result against CBSAs, commuting zones, and a modifiable areal unit problem (MAUP) ladder.

The pipeline is written in Stata 18 with Python helpers for the graph and clustering steps. `_master.do` is the only entry point.

## Requirements

**Stata 18** (StataNow or later) with these community packages installed from SSC:

| Package | Commands used | Purpose |
|---|---|---|
| `gtools` | `gcollapse`, `gegen`, `gcontract` | fast aggregation of the OD rows |
| `geoplot` and `geoframe` (Jann) | `geoframe`, `geoplot` | maps in `gen/` and `_inc/` |
| dependencies of `geoplot` | `palettes`, `colrspace`, `moremata` | pulled in by `geoplot` |

`spshape2dta` ships with Stata's Sp suite. `_master.do` checks for all of the commands in the first column and stops with a list of any that are missing.

**Python 3.9 or later** (developed on 3.14) with the packages in `_py/requirements.txt`: `pandas`, `numpy`, `scipy`, `scikit-learn`, `python-igraph`, `leidenalg`, `matplotlib`, and, for the contiguity step only, `geopandas` and `libpysal`.

```
python -m venv /path/to/python-env
/path/to/python-env/bin/pip install -r _py/requirements.txt
```

**Disk and network.** The LODES download is roughly 70 GB for both job-type universes, or 35 GB for the baseline universe alone. The processed intermediates add tens of gigabytes more. Two assessment scripts fetch American Community Survey (ACS) summary-file tables and the Census Gazetteer from `www2.census.gov` at run time and cache them beside their output, so the `assess` phase needs a network connection on first run.

**Shell.** Several scripts call Python as `shell VAR=1 "$py" script.py ...`, which is POSIX shell syntax. We have run the pipeline only on macOS and Linux. On Windows those calls would need to be rewritten.

## Inputs and directory layout

The pipeline reads from a project directory that sits outside this repository. We refer to it as `$proj`; its layout is fixed:

```
$proj/
  03_data/
    00_shapefiles/
      us_tract_2020/US_tract_2020.shp     2020 TIGER/Line tracts
      us_county_2020/US_county_2020.shp   2020 TIGER/Line counties
      us_state_2020/US_state_2020.shp     2020 TIGER/Line states
    01_raw/
      lodes/<st>_od_main_<JT>_<year>.csv.gz   LODES8 OD, main part
      lodes/<st>_od_aux_<JT>_<year>.csv.gz    LODES8 OD, aux part
      cbsa/list1_2023.xlsx                    OMB 2023 CBSA delineation
      tract_populationcenters_2020.csv        2020 Centers of Population, tract file
      czones/county_to_cz.dta (or .csv)       county-to-commuting-zone crosswalk (optional)
    02_processed/                             written by the pipeline
  05_output/                                  written by the pipeline
```

`grab/fetch_lodes.sh` downloads the LODES files and the CBSA workbook into `01_raw/` and prints the URLs for the three shapefile sets and the population-centers file, which are large enough that we leave them to a manual download:

```
export LEHD_PROJ=/path/to/project
bash grab/fetch_lodes.sh
```

Set `JOB_TYPES="JT00"` in the environment to fetch only the all-jobs universe. The primary-jobs universe (JT01) feeds one robustness table and nothing else. The population-centers file is the tract file `CenPop2020_Mean_TR.txt` from the Census Bureau, renamed to `tract_populationcenters_2020.csv`; county centroids are derived from it, so the county file is not needed. The commuting-zone crosswalk is not a Census product and is optional; `assess/assess_06_czones.do` skips its one table when the file is absent.

One reference file also sits in `grab/`. `chosen_params.csv` records the calibrated parameters the paper uses (tract and county), for readers who want them without running the calibration phase. The same values are set in `_inc/calib_params.do`, which is what the pipeline reads. Nothing in the code reads the copy in `grab/`.

## Configuration

`_master.do` takes its paths from environment variables and falls back to `/path/to/...` placeholders that must be edited if the variables are not set:

| Variable | Stata global | Meaning |
|---|---|---|
| `LEHD_PROJ` | `$proj` | project root described above |
| `LEHD_PROG` | `$program` | this repository (the folder containing `_master.do`) |
| `LEHD_OUT` | `$output` | output root; default `$proj/05_output` |
| `LEHD_PYENV` | `$pyenv` | Python virtual environment |
| `LEHD_PY` | `$py` | Python interpreter; default `$pyenv/bin/python` |

The remaining knobs are globals near the top of `_master.do`:

- `$yr_start`, `$yr_end`: pooling window, 2019 to 2023.
- `$levels`: geographies to cluster, `tract county`.
- `$job_types`: LODES universes to build, `JT00 JT01`.
- `$noncontig_fips`: state FIPS codes dropped before any computation (Alaska, Hawaii, and the territories).
- `$seed`: 202600331. Every stochastic step (Leiden, bootstrap stability, permutation nulls, random partitions) receives it.
- `$do_data` through `$do_export`, and the `$do_robust_*` switches: set a phase to 0 to skip it. The calibration sweeps and the validity checks have their own switches inside the `$do_calib` block.

Leiden keeps the best of N seeds per division. The calibration sweep sets N to 10; the production run in `cluster/cluster_02_leiden.do` uses the default of 5 unless `LEIDEN_NSEEDS` is exported. The calibration grids (alpha, distance cap, HC k, Leiden resolution) are defined in `_master.do` beside the calls that use them.

## Running

```
export LEHD_PROJ=/path/to/project
export LEHD_PROG=/path/to/lehd-commuting-clusters
export LEHD_PYENV=/path/to/python-env
stata-mp -b do _master.do
```

The calibration and robustness phases account for most of the run time. `_master.do` writes per-script run times to `$output/logs/run_times_<run_tag>_<date>.txt`, and every script opens its own log in `$output/logs/`.

## Pipeline

The phases run in the order below. Within a phase, the file number is the run order except where noted. Robustness reruns re-enter the cluster phase through `_inc/run_cluster.do` under a different `$run_tag`, so every intermediate is namespaced by that tag.

### data

| Script | Purpose |
|---|---|
| `data_01_state_lookup.do` | Build the 49-row state FIPS to abbreviation lookup used to enumerate LODES files. |
| `data_02_create_files.do` | Read the main and aux OD files for one job type and year; keep home block, work block, and total jobs. |
| `data_03_pool_multi.do` | Aggregate blocks to tract and county and pool the five years. |
| `data_04_cbsa_crosswalk.do` | Build the county and tract to CBSA crosswalk from the OMB workbook. |

### build

| Script | Purpose |
|---|---|
| `build_01_pflows.do` | Convert pooled OD counts into place-of-residence flow shares with resident labor force. Runs once per job type and once per minimum-flow cut (2, 3, 5, 10, 20) for the threshold rerun. |
| `build_02_centroids.do` | Population-weighted tract centroids from the Centers of Population file; county centroids derived from them. |
| `build_03_divisions.do` | Snake-scan the country into overlapping divisions so that HC runs on tractable blocks. Runs once per division geometry in `_inc/div_specs.do`; `baseline` is the one the paper uses. |

### calib

| Script | Purpose |
|---|---|
| `calib_01_alpha_cap.do` | Sweep the flow-versus-distance weight alpha and the distance cap. |
| `calib_03_hc_k.do` | Sweep the national HC cluster count k. Runs before `calib_02`; the two are independent. |
| `calib_02_resolution.do` | Sweep the Leiden CPM resolution gamma. |
| `calib_04_pick_params.do` | Persist the chosen alpha, cap, k, and gamma to `$temp/chosen_params_<lvl>.dta`. We chose the values by reading the sweep diagnostics and set them in `_inc/calib_params.do`; the cluster and robust phases read only these files. |
| `calib_05_validity.do` | Bootstrap stability (Hennig 2007) and geography-only and permutation-null benchmarks for the chosen cut. |

### cluster

| Script | Purpose |
|---|---|
| `cluster_01_hc.do` | Write the per-division dissimilarity matrices and run HC to the chosen k. |
| `cluster_02_leiden.do` | Run Leiden CPM at the chosen gamma on the same matrices. |
| `cluster_03_compile.do` | Align the HC and Leiden labelings on a common geoid frame, one row per unit per division. |
| `cluster_04_fuse.do` | Stitch division-local clusters into national ones. Two clusters in adjacent divisions are unioned when they share at least `min_shared` overlap units and that overlap is at least `frac_min` of the smaller cluster. |
| `cluster_05_reconcile.do` | Assign each unit that falls in the overlap of two divisions to one division by its flow score. |
| `cluster_06_harmonize.do` | Replace division-local cluster identifiers with the globally unique ones from the fusion maps. |
| `cluster_07_contiguity.do` | Enforce spatial contiguity using a queen adjacency built from the TIGER shapefiles. |
| `cluster_08_absorb_singletons.do` | Absorb one-unit clusters into their dominant commuting partner and write `$temp/final_assignments_<run_tag>_<lvl>.dta`. |

### robust

Each script sets a `$run_tag`, reruns the cluster phase, and restores the baseline namespace through `_inc/robust_prologue.do` and `_inc/robust_epilogue.do`.

| Script | Rerun |
|---|---|
| `robust_01_window_3yr.do` | 2021 to 2023 window. |
| `robust_02_window_1yr.do` | Each single year. |
| `robust_03_window_baseline_ref.do` | Placeholder for the five-year window, which is the baseline run itself; it reruns nothing. |
| `robust_05_threshold.do` | Minimum-flow cuts of 2, 3, 5, 10, and 20 workers. |
| `robust_06_params.do` | One-at-a-time perturbations of alpha, cap, k, and gamma from `_inc/sens_grid.do`. |
| `robust_07_divisions.do` | Alternative snake-scan geometries. |
| `robust_08_alpha_full.do` | Full pipeline at alpha values from `$alpha_full_grid` (default 0.0 and 0.2). |
| (in `_master.do`) | Primary jobs (JT01) universe. |

There is no `robust_04`; the numbering has a gap.

### assess

| Script | Purpose |
|---|---|
| `assess_01_elbow.do` | Elbow diagnostics for HC and Leiden. |
| `assess_02_agreement.do` | Agreement between the HC and Leiden partitions. |
| `assess_03_temporal.do` | Stability across the window reruns. |
| `assess_04_sensitivity.do` | Sensitivity across the parameter reruns. |
| `assess_05_cbsa.do` | Preservation of 2023 CBSAs. |
| `assess_06_czones.do` | Comparison with USDA/ERS commuting zones (optional input). |
| `assess_07_partition_quality.do`, `assess_08_partition_quality_alpha.do` | Silhouette and related quality measures, overall and by alpha. |
| `assess_09_maup.do`, `assess_10_maup_paramsweep.do` | MAUP ladder: a tract-level ACS relationship re-estimated at ten areal-unit definitions. |
| `assess_11_robinson.do` | Robinson-style ecological correlation across the same ladder; fetches the ACS pairs file. |
| `assess_12_random_null.do` | Random contiguous and non-contiguous partitions as a null for the ladder. Requires `assess_11` to have run. |

### gen, report, export

`gen/tab_*.do` and `gen/fig_*.do` write the manuscript tables (`.tex` and `.csv`) and figures (`.png`) to `$output/tables` and `$output/figures`; `report/report_01`, `_02`, and `_04` write the descriptive table, the worked dissimilarity example, and the validation table (there is no `report_03`). `export/report_10_export.do` writes the deliverables:

```
$output/tract_clusters_<run_tag>.dta / .csv
$output/county_clusters_<run_tag>.dta / .csv
```

Each row is a geoid with its final cluster, the HC and Leiden cluster identifiers that produced it, the snake-scan division, state FIPS, centroid coordinates, and population.

## Helpers

`_inc/` holds included fragments: the namespace defaults (`baseline_namespaces.do`), the timing wrapper (`timing_init.do`, `_timeit.do`, `timing_done.do`), the division specifications, the LaTeX table writer (`booktabs.do`), the map engines, and the robustness prologue and epilogue. `_py/` holds the Python entry points; each module docstring states the script that calls it and the arguments it expects.

## Citation

See `CITATION.cff`. The code is released under the University of Illinois/NCSA Open Source License (`LICENSE`).
