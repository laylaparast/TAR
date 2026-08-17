########################################################################
# tar_run_setting_081726.R
#
# Runs calibration (population truth) and the main simulation loop for
# ONE setting (1-8). Set `setting_to_run` before sourcing this file --
# this lets the eight settings run as eight independent jobs in
# parallel, e.g. via the eight run_setting<N>.R launcher scripts:
#
#   setting_to_run <- 1
#   source("tar_run_setting_081726.R")
#
# After all eight settings have finished and written their
# results_settingN_v3.rds files, run tar_combine_results_081726.R once
# to build the figures and summary tables.
########################################################################

if (!exists("setting_to_run") || !(setting_to_run %in% 1:8)) {
  stop("tar_run_setting_081726.R: `setting_to_run` must be set to an integer 1-8 before sourcing this file, e.g.:\n  setting_to_run <- 1\n  source(\"tar_run_setting_081726.R\")")
}

source("tar_functions_081726.R")

n0 <- 500
n1 <- 500
epsilon <- 0.05

###########################################################################
# Calibration: population truth for R_full and R_S0 for this setting.
#   - Settings 1,2,3,4,8 (linear, Normal candidates): linear_truth(), an
#     exact closed-form calculation.
#   - Settings 5,6,7 (nonlinear/heterogeneous): plugin_truth_nonlinear(),
#     a large-sample Monte Carlo plug-in calculation (n_big = 2e6 for
#     Setting 7, 1e6 for Settings 5-6).
###########################################################################

cat(sprintf("=== Calibration: population truth for Setting %d ===\n", setting_to_run))
s0_cols <- true_minimal_set[[as.character(setting_to_run)]]
if (setting_to_run %in% c(1, 2, 3, 4, 8)) {
  tr <- linear_truth(setting_to_run)
  Rf <- tr$R_full; Rs <- tr$R_S0
} else {
  full_cols <- paste0("s", 1:9)
  n_big_s <- if (setting_to_run == 7) 2e6 else 1e6
  Rf <- plugin_truth_nonlinear(setting_to_run, full_cols, n_big = n_big_s, seed = 100 + setting_to_run)$R
  Rs <- plugin_truth_nonlinear(setting_to_run, s0_cols,   n_big = n_big_s, seed = 200 + setting_to_run)$R
}
calib_row <- data.frame(
  Setting = setting_to_run, R_full = round(Rf, 4), R_S0 = round(Rs, 4),
  gap_S0_minus_full = round(Rs - Rf, 4), S0 = paste(s0_cols, collapse = ",")
)
print(calib_row)
if (abs(calib_row$gap_S0_minus_full) > epsilon) {
  warning(sprintf("Calibration gap exceeds epsilon (%.2f) for Setting %d -- R_S0 is not a trustworthy near-optimal reference for this setting.",
                   epsilon, setting_to_run))
}
cat("\n")

###########################################################################
# Main simulation loop for this setting.
###########################################################################

n.sims <- 500
setting_labels <- c(
  "1" = "Setting 1: baseline (rho=0)",
  "2" = "Setting 2: correlated (rho=0.5)",
  "3" = "Setting 3: rogue variables",
  "4" = "Setting 4: redundant / minimal set",
  "5" = "Setting 5: nonlinear baseline+corr",
  "6" = "Setting 6: nonlinear rogue",
  "7" = "Setting 7: nonlinear redundant/minimal set",
  "8" = "Setting 8: many-candidate (k=80)"
)

# Bootstrap: run.bootstrap = TRUE reruns the full estimation algorithm
# (including re-selection of lambda*) boot.num times per replicate, per
# approach. Recommended to sanity-check with small numbers first
# (run.bootstrap <- TRUE; n.sims <- 5; boot.num <- 20) before a full run.
run.bootstrap <- TRUE
boot.num <- 200

# Parallel backend for .bootstrap_var()'s boot.num resamples.
# parallelly::availableCores() (if installed) respects a cluster
# scheduler's actual core allocation rather than the whole node's
# hardware count; max_workers is a hard ceiling on top of that.
max_workers <- 32
n_workers <- if (requireNamespace("parallelly", quietly = TRUE)) {
  min(parallelly::availableCores(omit = 1), max_workers)
} else {
  min(max(1, parallel::detectCores() - 1), max_workers)
}
if (future::supportsMulticore()) {
  future::plan(future::multicore, workers = n_workers)
} else {
  future::plan(future::multisession, workers = n_workers)
}
cat(sprintf("Parallel backend: %s with %d workers.\n", class(future::plan())[1], n_workers))

setting <- setting_to_run
k <- if (setting == 8) 80 else 9
surrogate_names <- paste0("s", 1:k)

cat("Running", setting_labels[as.character(setting)], "...\n")
set.seed(1000 + setting)

R_full     <- numeric(n.sims)
R_smoothed <- numeric(n.sims)
R_tar      <- numeric(n.sims)
R_taronly  <- numeric(n.sims)

sel_smoothed <- matrix(0, n.sims, k, dimnames = list(NULL, surrogate_names))
sel_tar      <- matrix(0, n.sims, k, dimnames = list(NULL, surrogate_names))
sel_taronly  <- matrix(0, n.sims, k, dimnames = list(NULL, surrogate_names))

# Bootstrap SE (feeds the table's ASE column); NA throughout if
# run.bootstrap == FALSE.
R_full_se     <- rep(NA_real_, n.sims)
R_smoothed_se <- rep(NA_real_, n.sims)
R_tar_se      <- rep(NA_real_, n.sims)
R_taronly_se  <- rep(NA_real_, n.sims)

n_full_failures <- 0

for (jj in 1:n.sims) {
  dat <- simulate_setting(n1, n0, setting)

  full_r <- full.surrogate.approach(dat, surrogate_names = surrogate_names, var = run.bootstrap, boot.num = boot.num)
  R_full[jj] <- full_r$pte
  if (is.na(full_r$pte)) n_full_failures <- n_full_failures + 1

  rs <- smoothed.regularized(dat, surrogate_names = surrogate_names, nlambda = 20, var = run.bootstrap, boot.num = boot.num)
  R_smoothed[jj] <- rs$pte
  sel_smoothed[jj, names(rs$selection)] <- 1

  rtar <- tar.approach(dat, surrogate_names = surrogate_names, folds = 5, epsilon = epsilon, nlambda = 20, var = run.bootstrap, boot.num = boot.num)
  R_tar[jj] <- rtar$pte
  sel_tar[jj, names(rtar$selection)] <- 1

  rto <- tar.weightsonly(dat, surrogate_names = surrogate_names, nlambda = 20, var = run.bootstrap, boot.num = boot.num)
  R_taronly[jj] <- rto$pte
  sel_taronly[jj, names(rto$selection)] <- 1

  if (run.bootstrap) {
    R_full_se[jj]     <- full_r$var$R.s_se
    R_smoothed_se[jj] <- rs$var$R.s_se
    R_tar_se[jj]      <- rtar$var$R.s_se
    R_taronly_se[jj]  <- rto$var$R.s_se
  }

  if (jj %% 100 == 0) cat("  rep", jj, "/", n.sims, "\n")
}

if (n_full_failures > 0) {
  cat(sprintf("NOTE: full-surrogate approach failed on %d of %d replicates in Setting %d (k=%d, n1=%d).\n",
              n_full_failures, n.sims, setting, k, n1))
}

result_this_setting <- list(
  R_full = R_full, R_smoothed = R_smoothed, R_tar = R_tar, R_taronly = R_taronly,
  sel_smoothed = sel_smoothed, sel_tar = sel_tar, sel_taronly = sel_taronly,
  R_full_se = R_full_se, R_smoothed_se = R_smoothed_se, R_tar_se = R_tar_se, R_taronly_se = R_taronly_se,
  n_full_failures = n_full_failures
)

saveRDS(result_this_setting, paste0("results_setting", setting, "_v3.rds"))
cat(sprintf("\nSetting %d done. Saved results_setting%d_v3.rds.\n", setting, setting))
cat("Once all 8 settings have finished, run tar_combine_results_081726.R once to build the plots and tables.\n")
