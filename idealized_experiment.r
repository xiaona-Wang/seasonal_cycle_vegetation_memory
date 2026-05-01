## Idealized Experiment

generate_synthetic_signal <- function(
  Ntim,
  NpY,
  id,
  varphi,
  amp_mod = 0.30,
  phase_mod = pi / 12,
  baseline = 0.10,
  peak = 0.70,
  residual_sd = 0.05,
  h1 = 1.00,
  h2 = 0.30,
  h3 = 0.18,
  white_noise_sd = 0,
  trend_scale = 0.20
) {
  set.seed(id)
  if (!is.numeric(Ntim) || length(Ntim) != 1 || Ntim <= 0) {
    stop("")
  }
  if (!is.numeric(NpY) || length(NpY) != 1 || NpY <= 1) {
    stop("")
  }
  if (baseline >= peak) {
    stop("")
  }
  if (Ntim %% NpY != 0) {
    warning("")
  }
  n_year <- Ntim / NpY
  tt <- 0:(Ntim - 1)
  t_year <- tt / NpY
  x <- seq(0, 1, length.out = Ntim)
  z1 <- 1 / (1 + exp(-18 * (x - 0.18)))
  z2 <- 1 / (1 + exp(-22 * (x - 0.50)))
  z3 <- 1 / (1 + exp(-16 * (x - 0.80)))
  trend_raw <- 0.14 * z1 - 0.08 * z2 + 0.12 * z3 +
	0.02 * sin(2 * pi * 0.75 * x) +
	0.01 * cos(2 * pi * 0.35 * x)
  trend_raw <- trend_raw - min(trend_raw)
  trend <- trend_raw / max(trend_raw) * trend_scale
  amp_t <- 1 +
    amp_mod * sin(2 * pi * t_year / max(2, n_year / 1.8)) +
    0.12 * cos(2 * pi * t_year / max(2, n_year / 3.0))
  phase_t <- phase_mod * sin(2 * pi * t_year / max(2, n_year / 2.4)) +
    0.20 * phase_mod * cos(2 * pi * t_year / max(2, n_year / 4.2))
  seasonal_raw <- amp_t * (
    h1 * cos(2 * pi * 1 * t_year + phase_t - pi / 3) +
    h2 * cos(2 * pi * 2 * t_year + 0.6 * phase_t - 0.5) +
    h3 * cos(2 * pi * 3 * t_year + 0.3 * phase_t)
  )
  raw_min <- min(seasonal_raw)
  raw_max <- max(seasonal_raw)
  if (raw_max == raw_min) {
    sin_component <- rep((baseline + peak) / 2, Ntim)
  } else {
    sin_component <- baseline +
      (seasonal_raw - raw_min) / (raw_max - raw_min) * (peak - baseline)
  }
  make_ar1 <- function(n, phi, sd = 0.05) {
    e <- stats::rnorm(n, mean = 0, sd = 1)
    out <- numeric(n)
    out[1] <- e[1]
    for (i in 2:n) {
      out[i] <- phi * out[i - 1] + sqrt(1 - phi^2) * e[i]
    }
    out <- out / stats::sd(out) * sd
    out
  }
  residual <- make_ar1(Ntim, varphi, sd = residual_sd)
  X <- sin_component + trend + residual + stats::rnorm(Ntim, mean = 0, sd = white_noise_sd)
  return(list(
    X = X,
    sin_component = sin_component,
    trend = trend,
    residual = residual,
    components = list(
      seasonal_raw = seasonal_raw,
      amp_t = amp_t,
      phase_t = phase_t
    ),
    params = list(
      Ntim = Ntim,
      NpY = NpY,
      id = id,
      varphi = varphi,
      amp_mod = amp_mod,
      phase_mod = phase_mod,
      baseline = baseline,
      peak = peak,
      residual_sd = residual_sd,
      h1 = h1,
      h2 = h2,
      h3 = h3,
      trend_scale = trend_scale
    )
  ))
}

evaluate_synthetic_tau_test <- function(
  Ntim,
  NpY,
  id,
  varphi,
  L1 = 368,
  L2 = 69,
  t.window = 367,
  s.window = 47,
  trend_scale = 0.20,
  amp_mod = 0.30,
  phase_mod = pi / 12,
  white_noise_sd = 0,
  path
) {
  y <- generate_synthetic_signal(Ntim, NpY, id, varphi,
	 amp_mod = amp_mod,phase_mod = phase_mod,white_noise_sd = white_noise_sd,trend_scale = trend_scale)
  dssa <- decomp_ssa(y$X, frequency=NpY, L1=L1, L2=L2)
  dstl <- decomp_stl(y$X, frequency=NpY, t.window=t.window, s.window=s.window)
  dtrad <- decomp_trad(y$X, frequency=NpY)
  mat <- cbind(dssa$residual, dstl$residual, dtrad$residual, y$residual) %>%
    apply(2, autocorrelation_calc, floor(Ntim / 2.5))
  tau_results <- apply(mat, 2, function(ac) {
    result <- estimate_tau_lsq(ac, return_quality = TRUE)
    c(tau = result$tau, R2 = result$R2, RMSE = result$RMSE)
  })
  tau_df <- as.data.frame(t(tau_results))
  colnames(tau_df) <- c("tau", "R2", "RMSE")
  tau_df$tau[tau_df$R2 < 0.4 | tau_df$RMSE > 0.5 | !is.finite(tau_df$tau)] <- NA_real_
  tau_final <- tau_df$tau
  utils::write.csv(tau_final, path, row.names = FALSE)
  return(tau_final)
}
