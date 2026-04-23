source("step0_function.r")
library(parallel)
library(doParallel)
library(foreach)
library(dplyr)
path3 <- "/Volumes/STUDY1/topic1/output/step3/"
dir.create(file.path(path3, "all_metrics_each_run"), showWarnings = FALSE)
cl <- makeCluster(10)
registerDoParallel(cl)
dat_all_metrics <- foreach(
  x = 1:1000,
  .combine = "rbind",
  .export = c(
    "calc_rmse", "calc_acf_rmse", "residual_accuracy_tau",
    "evaluate_synthetic_all_metrics",
    "generate_synthetic_signal",
    "decomp_ssa", "decomp_stl", "decomp_trad",
    "fill_gap", "fill_gap_auto",
    "autocorrelation_calc", "estimate_tau_lsq"
  ),
  .packages = c("zoo", "dplyr", "doParallel", "Rssa", "forecast")
) %dopar% {
  evaluate_synthetic_all_metrics(
    Ntim          = 46 * 23,
    NpY           = 46,
    id            = x,
    varphi        = 0.6,
    L1        = 368,
    L2        = 69,
    t.window  = 367,
    s.window  = 47,
    path          = file.path(path3, "all_metrics_each_run", sprintf("%05d_res.csv", x))
  )
}
stopCluster(cl)
write.csv(dat_all_metrics, file.path(path3, "all_metrics_1000runs.csv"), row.names = FALSE)

summary_mean <- dat_all_metrics |>
  dplyr::summarise(dplyr::across(where(is.numeric), ~ mean(.x, na.rm = TRUE)))
summary_sd <- dat_all_metrics |>
  dplyr::summarise(dplyr::across(where(is.numeric), ~ stats::sd(.x, na.rm = TRUE)))
summary_median <- dat_all_metrics |>
  dplyr::summarise(dplyr::across(where(is.numeric), ~ stats::median(.x, na.rm = TRUE)))
write.csv(summary_mean,   file.path(path3, "all_metrics_mean.csv"),   row.names = FALSE)
write.csv(summary_sd,     file.path(path3, "all_metrics_sd.csv"),     row.names = FALSE)
write.csv(summary_median, file.path(path3, "all_metrics_median.csv"), row.names = FALSE)