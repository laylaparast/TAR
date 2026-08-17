# run_setting1.R -- launcher for cluster submission.
# Sets which setting to run, then sources tar_run_setting_081726.R.
# Submit run_setting1.R through run_setting8.R together to run all 8
# settings in parallel. After all 8 have finished, run
# tar_combine_results_081726.R once to build the figures and tables.

setting_to_run <- 1
source("tar_run_setting_081726.R")
