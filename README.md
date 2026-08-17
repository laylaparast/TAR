# Reproducibility materials: Treatment-Aware Regularization (TAR)

This folder contains the R code needed to reproduce the simulation study results reported in the paper. All files below should be placed together in a single directory before running anything.

## Files

| File | Purpose |
|---|---|
| `tar_functions_081726.R` | All estimation methods (Full surrogate, Robust regularization, TAR weights-only, TAR), the eight data-generating settings, population-truth calculations, and plotting/table-writing functions. Sourced by the other scripts; never run directly. |
| `tar_run_setting_081726.R` | Runs the calibration and simulation loop for a single setting. Sourced by the eight `run_setting<N>.R` launchers below; not run directly. |
| `run_setting1.R` ... `run_setting8.R` | Launcher scripts, one per setting. Each sets which setting to run and sources `tar_run_setting_081726.R`. |
| `tar_combine_results_081726.R` | Run once, after all eight settings have finished. Builds the manuscript's figures, summary tables, and the large-sample R_Q diagnostic. |

## Requirements

R (version 4.0 or later recommended) with the following packages installed:

```r
install.packages(c("glmnet", "dplyr", "tidyr", "ggplot2", "patchwork",
                    "Rsurrogate", "MASS", "future.apply", "future"))
```

`parallelly` is optional but recommended (used to detect available cores correctly on shared/cluster machines; falls back to `parallel::detectCores()` if not installed).

## How to reproduce the results

**Step 1: Run all eight settings.**

Each setting is independent and can be run in any order, or in parallel as separate jobs (e.g., separate cluster submissions). For each setting `N` from 1 to 8, open R in the directory containing these files and run:

```r
source("run_settingN.R")
```

(substituting the actual number for `N`, e.g. `source("run_setting1.R")`). Each script simulates 500 datasets, applies all four estimation approaches, and saves the results to `results_settingN_v3.rds`.

By default, each setting also runs a nonparametric bootstrap (200 resamples per replicate) to compute standard errors, which is the most time-consuming part of the simulation. To skip the bootstrap for a faster initial check, edit `run.bootstrap <- TRUE` to `run.bootstrap <- FALSE` near the top of `tar_run_setting_081726.R` before sourcing.

All eight settings must finish (i.e., all eight `results_settingN_v3.rds` files must exist in the directory) before moving to Step 2.

**Step 2: Combine results.**

Once all eight `results_settingN_v3.rds` files exist, run:

```r
source("tar_combine_results_081726.R")
```

This produces:

- `table_settings1to4.tex`, `table_settings5to8.tex` — LaTeX tables of PTE estimates and selection summaries (Tables in the manuscript)
- `figure_selection_1to4.pdf`, `figure_selection_5to8.pdf` — combined selection-probability figures (Figures in the manuscript)
- `table_setting8fp.tex`, `setting8_false_positive_counts.csv` — Setting 8 noise-candidate selection summary
- `table_RQ_diagnostic.csv` — large-sample R_Q and surrogate-index overlap diagnostic for Settings 5-7
- `tar_simulation_summary.csv` — the full summary table underlying both LaTeX tables

Note that the R_Q diagnostic step is the slowest part of this script (a large-sample kernel smoothing calculation); expect it to take several minutes.

## Reproducibility notes

- All random data generation and resampling is seeded within the scripts, so results should be identical across runs on the same machine/R version.
- The simulation uses `n1 = n0 = 500` (500 treated, 500 control) per replicate, 500 replicates per setting, and 200 bootstrap resamples per replicate, matching the manuscript.
- Running all eight settings with the bootstrap enabled is computationally intensive; using a multi-core machine is strongly recommended (the scripts parallelize the bootstrap automatically via the `future` package).
