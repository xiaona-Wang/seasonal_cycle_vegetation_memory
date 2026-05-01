## Parameter Selection

find_best_ssa_params <- function(
  Ntim, NpY, id, varphi,
  path = NULL,
  seas_cor_min = 0.90,
  trend_cor_min = 0.90,
  seas_rmse_q = 0.30,
  trend_rmse_q = 0.50,
  acf_tol = 0.01,
  tau_tol = 0.01,
  return_full = FALSE,
  return_thresholds = FALSE
) {
  sig <- generate_synthetic_signal(Ntim, NpY, id, varphi)
  L1_grid <- ceiling(seq(NpY * 2, Ntim / 2.6, NpY))
  L2_grid <- ceiling(seq(NpY, NpY * 4, NpY / 2))
  out <- list()
  k <- 1
  for (L1 in L1_grid) {
    for (L2 in L2_grid) {
      if (L1 <= L2) next
      fit <- tryCatch(
        decomp_ssa(sig$X, frequency = NpY, L1 = L1, L2 = L2),
        error = function(e) NULL
      )
      if (is.null(fit)) next
      out[[k]] <- data.frame(
        method     = "SSA",
        L1         = L1,
        L2         = L2,
        trend_cor  = suppressWarnings(cor(fit$trend,    sig$trend,         use = "pairwise.complete.obs")),
        trend_rmse = calc_rmse(sig$trend, fit$trend),
        seas_cor   = suppressWarnings(cor(fit$season,   sig$sin_component, use = "pairwise.complete.obs")),
        seas_rmse  = calc_rmse(sig$sin_component, fit$season),
        resid_cor  = suppressWarnings(cor(fit$residual, sig$residual,      use = "pairwise.complete.obs")),
        resid_rmse = calc_rmse(sig$residual, fit$residual),
        delta_tau  = residual_accuracy_tau(sig$residual, fit$residual),
        acf_rmse   = calc_acf_rmse(sig$residual, fit$residual)
      )
      k <- k + 1
    }
  }
  res <- dplyr::bind_rows(out)
  if (!is.null(path)) {
    utils::write.csv(res, path, row.names = FALSE)
  }
  required_cols <- c(
    "seas_cor", "seas_rmse",
    "trend_cor", "trend_rmse",
    "acf_rmse", "delta_tau",
    "resid_cor", "resid_rmse"
  )
  miss_cols <- setdiff(required_cols, names(res))
  if (length(miss_cols) > 0) {
    stop("", paste(miss_cols, collapse = ", "))
  }
  dat <- res |>
    dplyr::filter(
      is.finite(seas_cor),
      is.finite(seas_rmse),
      is.finite(trend_cor),
      is.finite(trend_rmse),
      is.finite(acf_rmse),
      is.finite(delta_tau),
      is.finite(resid_cor),
      is.finite(resid_rmse)
    )
  if (nrow(dat) == 0) return(NULL)
  seas_rmse_max <- as.numeric(stats::quantile(dat$seas_rmse, probs = seas_rmse_q, na.rm = TRUE))
  trend_rmse_max <- as.numeric(stats::quantile(dat$trend_rmse, probs = trend_rmse_q, na.rm = TRUE))
  dat1 <- dat |>
    dplyr::filter(
      seas_cor >= seas_cor_min,
      seas_rmse <= seas_rmse_max
    )
  if (nrow(dat1) == 0) return(NULL)
  dat2 <- dat1 |>
    dplyr::filter(
      trend_cor >= trend_cor_min,
      trend_rmse <= trend_rmse_max
    )
  if (nrow(dat2) == 0) return(NULL)
  best_acf <- min(dat2$acf_rmse, na.rm = TRUE)
  dat3 <- dat2 |>
    dplyr::filter(acf_rmse <= best_acf + acf_tol)
  if (nrow(dat3) == 0) return(NULL)
  best_tau <- min(abs(dat3$delta_tau), na.rm = TRUE)
  dat4 <- dat3 |>
  dplyr::filter(abs(delta_tau) <= best_tau + tau_tol)
  if (nrow(dat4) == 0) return(NULL)
  selected <- dat4 |>
    dplyr::arrange(
      dplyr::desc(resid_cor),
      resid_rmse,
      acf_rmse,
      abs(delta_tau)
    )
  out_list <- list(best = selected[1, ])
  if (return_full) {
    out_list$all <- res
    out_list$selected <- selected
  }
  if (return_thresholds) {
    out_list$thresholds <- list(
      seas_cor_min   = seas_cor_min,
      trend_cor_min  = trend_cor_min,
      seas_rmse_q    = seas_rmse_q,
      trend_rmse_q   = trend_rmse_q,
      seas_rmse_max  = seas_rmse_max,
      trend_rmse_max = trend_rmse_max,
      acf_tol        = acf_tol,
      tau_tol        = tau_tol
    )
  }
  out_list
}

find_best_stl_params <- function(
  Ntim, NpY, id, varphi,
  path = NULL,
  seas_cor_min = 0.90,
  trend_cor_min = 0.90,
  seas_rmse_q = 0.30,
  trend_rmse_q = 0.50,
  acf_tol = 0.01,
  tau_tol = 0.01,
  return_full = FALSE,
  return_thresholds = FALSE
) {
  sig <- generate_synthetic_signal(Ntim, NpY, id, varphi)
  t.window.grid <- ceiling(seq(NpY, NpY * 4, NpY / 2)) * 2 - 1
  s.window.grid <- ceiling(seq(NpY, NpY * 4, NpY / 2)) + 1
  out <- list()
  k <- 1
  for (t.window in t.window.grid) {
    for (s.window in s.window.grid) {
      if (t.window <= s.window) next
      fit <- tryCatch(
        decomp_stl(sig$X, frequency = NpY, t.window = t.window, s.window = s.window),
        error = function(e) NULL
      )
      if (is.null(fit)) next
      out[[k]] <- data.frame(
        method     = "STL",
        t.window   = t.window,
        s.window   = s.window,
        trend_cor  = suppressWarnings(cor(fit$trend,    sig$trend,         use = "pairwise.complete.obs")),
        trend_rmse = calc_rmse(sig$trend, fit$trend),
        seas_cor   = suppressWarnings(cor(fit$season,   sig$sin_component, use = "pairwise.complete.obs")),
        seas_rmse  = calc_rmse(sig$sin_component, fit$season),
        resid_cor  = suppressWarnings(cor(fit$residual, sig$residual,      use = "pairwise.complete.obs")),
        resid_rmse = calc_rmse(sig$residual, fit$residual),
        delta_tau  = residual_accuracy_tau(sig$residual, fit$residual),
        acf_rmse   = calc_acf_rmse(sig$residual, fit$residual)
      )
      k <- k + 1
    }
  }
  res <- dplyr::bind_rows(out)
  if (!is.null(path)) {
    utils::write.csv(res, path, row.names = FALSE)
  }
  required_cols <- c(
    "seas_cor", "seas_rmse",
    "trend_cor", "trend_rmse",
    "acf_rmse", "delta_tau",
    "resid_cor", "resid_rmse"
  )
  miss_cols <- setdiff(required_cols, names(res))
  if (length(miss_cols) > 0) {
    stop("", paste(miss_cols, collapse = ", "))
  }
  dat <- res |>
    dplyr::filter(
      is.finite(seas_cor),
      is.finite(seas_rmse),
      is.finite(trend_cor),
      is.finite(trend_rmse),
      is.finite(acf_rmse),
      is.finite(delta_tau),
      is.finite(resid_cor),
      is.finite(resid_rmse)
    )
  if (nrow(dat) == 0) return(NULL)
  seas_rmse_max <- as.numeric(stats::quantile(dat$seas_rmse, probs = seas_rmse_q, na.rm = TRUE))
  trend_rmse_max <- as.numeric(stats::quantile(dat$trend_rmse, probs = trend_rmse_q, na.rm = TRUE))
  dat1 <- dat |>
    dplyr::filter(
      seas_cor >= seas_cor_min,
      seas_rmse <= seas_rmse_max
    )
  if (nrow(dat1) == 0) return(NULL)
  dat2 <- dat1 |>
    dplyr::filter(
      trend_cor >= trend_cor_min,
      trend_rmse <= trend_rmse_max
    )
  if (nrow(dat2) == 0) return(NULL)
  best_acf <- min(dat2$acf_rmse, na.rm = TRUE)
  dat3 <- dat2 |>
    dplyr::filter(acf_rmse <= best_acf + acf_tol)
  if (nrow(dat3) == 0) return(NULL)
  best_tau <- min(abs(dat3$delta_tau), na.rm = TRUE)
  dat4 <- dat3 |>
  dplyr::filter(abs(delta_tau) <= best_tau + tau_tol)
  if (nrow(dat4) == 0) return(NULL)
  selected <- dat4 |>
    dplyr::arrange(
      dplyr::desc(resid_cor),
      resid_rmse,
      acf_rmse,
      abs(delta_tau)
    )
  out_list <- list(best = selected[1, ])
  if (return_full) {
    out_list$all <- res
    out_list$selected <- selected
  }
  if (return_thresholds) {
    out_list$thresholds <- list(
      seas_cor_min   = seas_cor_min,
      trend_cor_min  = trend_cor_min,
      seas_rmse_q    = seas_rmse_q,
      trend_rmse_q   = trend_rmse_q,
      seas_rmse_max  = seas_rmse_max,
      trend_rmse_max = trend_rmse_max,
      acf_tol        = acf_tol,
      tau_tol        = tau_tol
    )
  }
  out_list
}