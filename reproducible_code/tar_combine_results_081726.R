########################################################################
# tar_combine_results_081726.R
#
# Run once, after all 8 run_setting<N>.R jobs have finished and written
# their results_settingN_v3.rds files. Builds the calibration table, the
# selection-probability figures, the two manuscript summary tables, and
# the large-sample R_Q / Q-overlap diagnostic (Settings 5-7).
#
# USAGE: source() or Rscript this file directly -- no setting_to_run
# needed.
########################################################################

source("tar_functions_081726.R")

epsilon <- 0.05

###########################################################################
# Calibration: population truth for R_full and R_S0, all 8 settings.
###########################################################################

cat("=== Calibration: population truth for R_full and R_{S0} ===\n")
calib <- data.frame()
for (s in 1:8) {
  s0_cols <- true_minimal_set[[as.character(s)]]
  if (s %in% c(1, 2, 3, 4, 8)) {
    tr <- linear_truth(s)
    Rf <- tr$R_full; Rs <- tr$R_S0
  } else {
    full_cols <- paste0("s", 1:9)
    n_big_s <- if (s == 7) 2e6 else 1e6
    Rf <- plugin_truth_nonlinear(s, full_cols, n_big = n_big_s, seed = 100 + s)$R
    Rs <- plugin_truth_nonlinear(s, s0_cols,   n_big = n_big_s, seed = 200 + s)$R
  }
  calib <- rbind(calib, data.frame(
    Setting = s, R_full = round(Rf, 4), R_S0 = round(Rs, 4),
    gap_S0_minus_full = round(Rs - Rf, 4), S0 = paste(s0_cols, collapse = ",")
  ))
}
print(calib)
if (any(abs(calib$gap_S0_minus_full) > epsilon)) {
  bad <- calib$Setting[abs(calib$gap_S0_minus_full) > epsilon]
  warning(sprintf("Calibration gap exceeds epsilon (%.2f) for Setting(s) %s.",
                   epsilon, paste(bad, collapse = ", ")))
}
cat("\n")

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

###########################################################################
# Load the 8 per-setting results.
###########################################################################

results_9 <- list()
for (setting in 1:8) {
  f <- paste0("results_setting", setting, "_v3.rds")
  if (!file.exists(f)) {
    stop(sprintf("%s not found -- all 8 results_settingN_v3.rds files must exist before running this script.", f))
  }
  results_9[[setting]] <- readRDS(f)
}
cat("Loaded all 8 results_settingN_v3.rds files.\n\n")

###########################################################################
# Figures: selection-probability panels for all 8 settings.
###########################################################################

panels <- vector("list", 7)

for (setting in 1:7) {
  res <- results_9[[setting]]
  props <- bind_rows(
    as.data.frame(res$sel_smoothed) %>% summarise(across(everything(), mean)) %>% mutate(Method = "Robust reg."),
    as.data.frame(res$sel_tar)      %>% summarise(across(everything(), mean)) %>% mutate(Method = "TAR"),
    as.data.frame(res$sel_taronly)  %>% summarise(across(everything(), mean)) %>% mutate(Method = "TAR (weights only)")
  ) %>%
    pivot_longer(cols = starts_with("s"), names_to = "Surrogate", values_to = "Selection_Prob") %>%
    mutate(Surrogate = factor(Surrogate, levels = paste0("s", 1:9)),
           Method = factor(Method, levels = c("Robust reg.", "TAR (weights only)", "TAR")))

  s0 <- true_minimal_set[[as.character(setting)]]

  p_standalone <- plot_selection_grid(props, setting_labels[as.character(setting)], s0_cols = s0)
  ggsave(paste0("selection_setting", setting, ".pdf"), p_standalone, width = 7, height = 3.2)

  panels[[setting]] <- plot_selection_grid(props, setting_labels[as.character(setting)], s0_cols = s0, compact = TRUE)
}

res8 <- results_9[[8]]
true3 <- c("s1", "s2", "s3")
noise_cols <- setdiff(colnames(res8$sel_smoothed), true3)

sel_true_df <- bind_rows(
  data.frame(Method = "Robust reg.",         Surrogate = true3, Selection_Prob = colMeans(res8$sel_smoothed[, true3, drop = FALSE])),
  data.frame(Method = "TAR (weights only)",  Surrogate = true3, Selection_Prob = colMeans(res8$sel_taronly[, true3, drop = FALSE])),
  data.frame(Method = "TAR",                 Surrogate = true3, Selection_Prob = colMeans(res8$sel_tar[, true3, drop = FALSE]))
)

sel_noise_long <- bind_rows(
  data.frame(Method = "Robust reg.",        Selection_Prob = colMeans(res8$sel_smoothed[, noise_cols])),
  data.frame(Method = "TAR (weights only)", Selection_Prob = colMeans(res8$sel_taronly[, noise_cols])),
  data.frame(Method = "TAR",                Selection_Prob = colMeans(res8$sel_tar[, noise_cols]))
)

p8_standalone <- plot_selection_setting8(sel_true_df, sel_noise_long)
ggsave("selection_setting8.pdf", p8_standalone, width = 8, height = 3.5)
panel8 <- plot_selection_setting8(sel_true_df, sel_noise_long, compact = TRUE)

fig_1to4 <- combined_selection_figure_1to4(panels[1:4])
ggsave("figure_selection_1to4.pdf", fig_1to4, width = 13, height = 8)

fig_5to8 <- combined_selection_figure_5to8(panels[5:7], panel8)
ggsave("figure_selection_5to8.pdf", fig_5to8, width = 13, height = 8)

cat("Saved figure_selection_1to4.pdf and figure_selection_5to8.pdf.\n")

fp_count <- write_setting8fp_table(res8$sel_smoothed, res8$sel_taronly, res8$sel_tar, noise_cols,
                                    file = "table_setting8fp.tex")
print(fp_count)
write.csv(fp_count, "setting8_false_positive_counts.csv", row.names = FALSE)

###########################################################################
# Summary tables: R_full, R_S0, Est, ASE, ESE, P(S0^C selected), mean #
# selected, for all 4 methods x 8 settings.
###########################################################################

summary_all <- data.frame()

for (setting in 1:8) {
  res <- results_9[[setting]]
  s0_cols <- true_minimal_set[[as.character(setting)]]
  R_S0   <- calib$R_S0[calib$Setting == setting]
  R_full <- calib$R_full[calib$Setting == setting]

  row_full     <- summarize_results(res$R_full,     NULL,             R_full, R_S0, s0_cols, "Full surrogate",     ase_vec = res$R_full_se)
  row_smoothed <- summarize_results(res$R_smoothed, res$sel_smoothed, R_full, R_S0, s0_cols, "Robust reg.",        ase_vec = res$R_smoothed_se)
  row_taronly  <- summarize_results(res$R_taronly,  res$sel_taronly,  R_full, R_S0, s0_cols, "TAR (weights only)", ase_vec = res$R_taronly_se)
  row_tar      <- summarize_results(res$R_tar,      res$sel_tar,      R_full, R_S0, s0_cols, "TAR",                ase_vec = res$R_tar_se)

  block <- rbind(row_full, row_smoothed, row_taronly, row_tar)
  block$Setting <- setting
  summary_all <- rbind(summary_all, block)
}

summary_all <- summary_all[, c("Setting", "Method", "R_full", "R_S0", "Est", "ASE", "ESE",
                                "P_S0C_selected", "Mean_num_selected")]
print(summary_all)
write.csv(summary_all, "tar_simulation_summary.csv", row.names = FALSE)

write_grouped_table(summary_all, settings = 1:4, setting_labels = setting_labels, file = "table_settings1to4.tex")
write_grouped_table(summary_all, settings = 5:8, setting_labels = setting_labels, file = "table_settings5to8.tex")

###########################################################################
# Large-sample R_Q / R_{Q0} diagnostic (Settings 5-7). Runs last: it is
# by far the slowest block here (a pairwise kernel smooth inside
# large_sample_RQ() -- see that function's header comment in
# tar_functions_081726.R), so running it after the figures and tables are
# already saved means it cannot block those from being produced.
###########################################################################

cat("=== Large-sample R_Q / R_{Q0} diagnostic (Settings 5-7) ===\n")
rq_diag <- data.frame()
for (s in 5:7) {
  full_cols <- paste0("s", 1:9)
  s0_cols   <- true_minimal_set[[as.character(s)]]

  rq_full <- large_sample_RQ(s, full_cols, seed = 300 + s)
  rq_s0   <- large_sample_RQ(s, s0_cols,   seed = 400 + s)

  rq_diag <- rbind(rq_diag, data.frame(
    Setting = s, Target = "R_Q (full, 9 candidates)",
    R_Q = round(rq_full$R_Q, 4), Delta = round(rq_full$Delta, 4),
    frac_outside_range_pop = round(rq_full$overlap_population$frac_outside_range, 5),
    frac_outside_1_99_pop  = round(rq_full$overlap_population$frac_outside_1_99, 5),
    overlap_finite_n500    = round(rq_full$overlap_finite_n500, 5)
  ))
  rq_diag <- rbind(rq_diag, data.frame(
    Setting = s, Target = "R_{Q0} (S0 only)",
    R_Q = round(rq_s0$R_Q, 4), Delta = round(rq_s0$Delta, 4),
    frac_outside_range_pop = round(rq_s0$overlap_population$frac_outside_range, 5),
    frac_outside_1_99_pop  = round(rq_s0$overlap_population$frac_outside_1_99, 5),
    overlap_finite_n500    = round(rq_s0$overlap_finite_n500, 5)
  ))
}
print(rq_diag)
write.csv(rq_diag, "table_RQ_diagnostic.csv", row.names = FALSE)
cat("Saved table_RQ_diagnostic.csv.\n\n")

cat("\nDone. Key outputs:\n")
cat("  - calib (printed above) : population truth check\n")
cat("  - table_RQ_diagnostic.csv : large-sample R_Q/R_{Q0} + Q-overlap diagnostic, Settings 5-7\n")
cat("  - selection_setting[1-8].pdf : standalone selection-probability panels\n")
cat("  - figure_selection_1to4.pdf, figure_selection_5to8.pdf : the two figures for the manuscript\n")
cat("  - tar_simulation_summary.csv : full summary table (R_full/R_S0/Est/ASE/ESE/P(S0^C selected)/Mean # selected)\n")
cat("  - table_settings1to4.tex, table_settings5to8.tex : LaTeX tabular fragments for the manuscript's two tables\n")
cat("  - table_setting8fp.tex, setting8_false_positive_counts.csv : mean # noise vars selected (of 77) in Setting 8\n")
