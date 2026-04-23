### functions used in the study
## =========================Preprocessing & Filtering===========================
library(Rssa)
library(zoo)
library(terra)
library(dplyr)
library(forecast)
fill_gap_auto <- function(x,
                          freq     = NULL,
                          L        = NULL,
                          groups   = 1:6,
                          fallback = c("interp", "none")) {

  fallback <- match.arg(fallback)
  x <- as.numeric(x)
  n <- length(x)
  if (n == 0 || all(is.na(x))) return(x)
  na_idx  <- is.na(x)
  n_na    <- sum(na_idx)
  n_valid <- n - n_na
  if (n_valid < 10) {
    warning("Too few non-NA values for SSA; using fallback.")
    if (fallback == "interp") {
      return(zoo::na.approx(x, na.rm = FALSE))
    } else {
      return(x)
    }
  }
  idx_non <- which(!is.na(x))
  x_interp <- stats::approx(
    x    = idx_non,
    y    = x[idx_non],
    xout = seq_len(n),
    rule = 2
  )$y
  if (is.null(L)) {
    if (!is.null(freq) && !is.na(freq)) {
      L_cand <- min(3 * freq, floor(n / 2))
    } else {
      L_cand <- floor(n / 2.5)
    }
    L <- max(10L, min(as.integer(L_cand), n - 2L))
  } else {
    L <- max(10L, min(as.integer(L), n - 2L))
  }
  rec <- tryCatch({
    ss  <- Rssa::ssa(x_interp, L = L, kind = "1d-ssa")
    Rssa::reconstruct(ss, groups = list(groups))[[1]]
  }, error = function(e) {
    message("SSA failed in fill_gap_auto: ", conditionMessage(e))
    NULL
  })
  if (is.null(rec)) {
    if (fallback == "interp") {
      return(zoo::na.approx(x, na.rm = FALSE))
    } else {
      return(x)
    }
  }
  z <- x
  z[na_idx] <- rec[na_idx]
  return(z)
}

fill_gap <- function(x, method) {
  x <- as.numeric(x)
  z <- rep(NA_real_, length(x))
  if (method == "interp") {
    z <- forecast::na.interp(x)
  } else if (method == "spline") {
    z <- zoo::na.spline(x)
  } else if (method == "ssa") {
    z <- fill_gap_auto(x)
  } else {
    z <- rep(NA_real_, length(x))
  }
  return(z)
}


flatten_2d_to_1d <- function(x) {
  stopifnot(length(dim(x)) == 2L)
  y <- numeric(dim(x)[1] * dim(x)[2])
  for (i in 1:dim(x)[1]) {
    y[((i - 1) * dim(x)[2] + 1):(i * dim(x)[2])] <- as.numeric(x[i, ])
  }
  return(y)
}

detect_vege_pixels <- function(
  tif,
  veg_codes      = 1:11,
  min_valid_frac = 0.75,
  cores          = 5,
  out            = NULL,
  wopt           = list(datatype="INT1U",
                        gdal=c("COMPRESS=LZW","NUM_THREADS=ALL_CPUS","BIGTIFF=YES")),
  .set_options   = TRUE,
  tempdir        = NULL,
  memfrac        = 0.8
) {
  stopifnot(inherits(tif, "SpatRaster"))
  if (.set_options) {
    terraOptions(memfrac = memfrac)
    if (!is.null(tempdir)) terraOptions(tempdir = tempdir)
  }
  f <- function(r) {
    if (all(is.na(r))) return(NA_real_)
    valid_frac <- sum(!is.na(r)) / length(r)
    if (valid_frac < min_valid_frac) return(NA_real_)
    vals <- r[!is.na(r)]
    if (any(!(vals %in% veg_codes))) return(0)
    if (length(vals) == 0) return(NA_real_)
    if (all(vals == vals[1])) return(1) else return(0)
  }
  app(tif, fun = f, cores = cores, filename = out, overwrite = TRUE, wopt = wopt)
}

detect_flat_pixels <- function(
  tif,
  min_valid_frac = 0.75,
  thr_mean_high  = 0.5,
  thr_mean_low   = 0.03,
  thr_sd         = 0.0001,
  thr_dx         = 0.0001,
  run_len        = 13 - 1,
  cores          = 1,
  out            = NULL,
  wopt           = list(datatype="INT1U", gdal=c("COMPRESS=LZW","NUM_THREADS=ALL_CPUS","BIGTIFF=YES")),
  .set_options   = TRUE,
  tempdir        = NULL,
  memfrac        = 0.8
) {
  stopifnot(inherits(tif, "SpatRaster"))
if (.set_options) {
  terraOptions(memfrac = memfrac)
  if (!is.null(tempdir)) {
    if (!dir.exists(tempdir)) dir.create(tempdir, recursive = TRUE)
    terraOptions(tempdir = tempdir)
  }
}
  f <- function(v) {
    n   <- if (is.null(dim(v))) 1L else nrow(v)
    out <- rep(NA_real_, n)
    keep <- rowSums(!is.na(v)) / ncol(v) >= min_valid_frac
    if (!any(keep)) return(out)
    vv <- v[keep, , drop = FALSE]
    m <- rowMeans(vv, na.rm = TRUE)
    s <- apply(vv, 1, sd, na.rm = TRUE)
    res <- rep(0, nrow(vv))
    bad <- is.na(m) | (m == 0) | is.na(s) | (s == 0)
    res[bad] <- NA
    res[!bad & (m > thr_mean_high)] <- 0
    res[!bad & (m < thr_mean_low)]  <- 1
    res[!bad & (s <= thr_sd)]       <- 1
    undecided <- which(!bad & !( (m < thr_mean_low) | (s <= thr_sd) | (m > thr_mean_high) ))
    if (length(undecided) > 0) {
      vv2 <- vv[undecided, , drop = FALSE]
      dx  <- t(apply(vv2, 1, diff))
      dx[is.na(dx)] <- thr_dx + 1
      stable <- abs(dx) < thr_dx
      long_run <- apply(stable, 1, function(z) { r <- rle(z); any(r$values & (r$lengths >= run_len)) })
      res[undecided[long_run]] <- 1
    }
    out[keep] <- res
    out
  }
  app(tif, fun = f, cores = cores, filename = out, overwrite = TRUE, wopt = wopt)
}

# ==============================Decomposition==================================
decomp_ssa <- function(x,
                       frequency,
					             L1,
                       L2) {
    if (length(stats::na.omit(x)) < (length(x) * 0.75)) {
        return(list(
		  trend = rep(NA_real_, length(x)),
		  season = rep(NA_real_, length(x)),
		  residual = rep(NA_real_, length(x))
		))
    }
	datv <- as.numeric(x)|>fill_gap('interp')
    meanvalue <- mean(datv, na.rm = TRUE)
    datv <- datv - meanvalue
    ss1 <- Rssa::ssa(datv, L = L1)
	g1_pgram <- grouping.auto(ss1, grouping.method = "pgram",base = "series",
	  freq.bins = list(
		trend = c(0, 1/(frequency*2)),
		other = c(1/(frequency*2), 0.5)
	  ),
	  threshold = 0,method = "constant"
	)
	r1_pgram <- reconstruct(ss1, groups = g1_pgram)
	trend1 <- if ("trend" %in% names(r1_pgram)) {
	  r1_pgram$trend
	} else {
	  reconstruct(ss1, groups = list(trend = 1:2))$trend
	}
	res1 <- datv - trend1
	ss2 <- Rssa::ssa(res1, L = L2)
	g2 <- grouping.auto(
	  ss2,
	  grouping.method = "pgram",
	  base = "series",
	freq.bins = list(
	  lowleft    = c(0,1/(frequency*2)),  
	  annual     = c(1/(frequency*2), 3/(frequency*2)),
	  semiannual = c(3/(frequency*2), 5/(frequency*2)),
	  triannual  = c(5/(frequency*2), 7/(frequency*2)),
	  highfreq   = c(7/(frequency*2), 0.5)
	),
	  threshold = 0, method = "constant"
	)
	r2 <- reconstruct(ss2, groups = g2)
	nz <- function(x) if (is.null(x)) 0 else x
	annual     <- nz(r2$annual)
	semiannual <- nz(r2$semiannual)
	triannual  <- nz(r2$triannual)
	seasonal <- annual + semiannual + triannual
	residual <- res1 - seasonal
	trend_final <- trend1 + meanvalue
    return(list(
        trend = as.numeric(trend_final) ,
        season = as.numeric(seasonal),
        residual = as.numeric(residual)
    ))
}

## mean seasonal cycle, y (sample number per year)
mean_cycle_calc <- function(x, y) {
  x <- as.numeric(x)
  n <- length(x)
  m <- ceiling(n / y)
  x_pad <- rep(NA_real_, y * m)
  x_pad[seq_len(n)] <- x
  dim(x_pad) <- c(y, m)
  mu <- rowMeans(x_pad, na.rm = TRUE)
  rep(mu, length.out = n)
}

decomp_trad <- function(x, frequency) {
  if (length(stats::na.omit(x)) < (length(x) * 0.75)) {
    return(data.frame(trend = NA_real_, season = NA_real_, residual = NA_real_))
  }
  datv <- as.numeric(x)|>fill_gap('interp')
  valid <- is.finite(datv)
  if (sum(valid) < 3) {
    return(data.frame(trend = NA_real_, season = NA_real_, residual = NA_real_))
  }
  tr_fit <- stats::lm(datv[valid] ~ seq_along(datv)[valid])
  trend_v <- rep(NA_real_, length(datv))
  trend_v[valid] <- stats::fitted(tr_fit)
  x_detrend <- datv - trend_v
  seas_v <- mean_cycle_calc(x_detrend, frequency)
  remainder <- x_detrend - seas_v
  data.frame(trend = trend_v , season = seas_v, residual = remainder)
}

## stl method
decomp_stl <- function(x, frequency, t.window, s.window) {
  if (length(stats::na.omit(x)) < (length(x) * 0.75)) {
    return(data.frame(trend = NA_real_, season = NA_real_, residual = NA_real_))
  }
  datv <- as.numeric(x)|>fill_gap('interp')
  datv[!is.finite(datv)] <- mean(datv, na.rm = TRUE)
  tsx <- stats::ts(datv, frequency = frequency)
  dat <- stats::stl(tsx, s.window = as.numeric(s.window), t.window = t.window, l.window = NULL)
  if (is.null(dat)) {
    return(data.frame(trend = NA_real_, season = NA_real_, residual = NA_real_))
  }
  data.frame(
    trend    = dat$time.series[, "trend"] ,
    season   = dat$time.series[, "seasonal"],
    residual = dat$time.series[, "remainder"]
  )
}

calc_rmse <- function(ref, test) {
  ok <- is.finite(ref) & is.finite(test)
  ref <- ref[ok]
  test <- test[ok]
  if (length(ref) == 0) return(NA_real_)
  sqrt(mean((test - ref)^2, na.rm = TRUE))
}

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

## ==============================Metrics & Evaluation===========================
# Compute autocorrelation up to lag l
autocorrelation_calc <- function(x, l) {
  x <- na.omit(as.numeric(x))
  n <- length(x)
  if (n < 2L) return(rep(NA_real_, l + 1))
  mean_x <- mean(x)
  x_dev <- x - mean_x
  d <- sum(x_dev^2)
  if (d == 0) return(rep(NA_real_, l + 1))
  res <- numeric(l + 1)
  for (t in 0:l) {
    res[t + 1] <- sum(x_dev[1:(n - t)] * x_dev[(1 + t):n]) / d * (n / (n - t))
  }
  res
}

## calculate tau
estimate_tau_lsq <- function(x, cutoff = 0.1, frac = 2.5, tau_max = 200, return_quality = TRUE) {
  x <- x[is.finite(x) & x > 0]
  n <- length(x)
  if (n < 4L) {
    return(if (return_quality) list(tau = NA_real_, R2 = NA_real_, RMSE = NA_real_) else NA_real_)
  }
  neg_pos   <- which(x < 0)[1]
  below_pos <- which(x < cutoff)[1]
  istop <- suppressWarnings(min(c(floor(n / frac), neg_pos, below_pos), na.rm = TRUE))
  if (!is.finite(istop) || istop < 3L) istop <- n
  data  <- x[1:istop]
  data  <- data[data > 0]
  if (length(unique(data)) < 3) {
    return(if (return_quality) list(tau = NA_real_, R2 = NA_real_, RMSE = NA_real_) else NA_real_)
  }
  tdata <- seq_len(length(data))
  log_y <- log(data)
  lmfit <- lm(log_y ~ tdata)
  slope <- coef(lmfit)[2]
  if (!is.finite(slope) || slope >= 0) {
    return(if (return_quality) list(tau = NA_real_, R2 = NA_real_, RMSE = NA_real_) else NA_real_)
  }
  tau_init <- -1 / slope
  A_init   <- exp(coef(lmfit)[1])
  tau_final <- NA_real_
  R2_final <- NA_real_
  RMSE_final <- NA_real_
  if (tau_init < tau_max && tau_init > 0.1) {
    fit <- tryCatch(
      nls(data ~ A * exp(-tdata / tau),
          start = list(A = A_init, tau = tau_init),
          control = nls.control(maxiter = 30, warnOnly = TRUE, minFactor = 1/2048)),
      error = function(e) NULL
    )
    if (!is.null(fit)) {
      tau <- as.numeric(coef(fit)["tau"])
      if (is.finite(tau) && tau > 0 && tau < tau_max) {
        tau_final <- tau
        pred <- predict(fit)
        ss_res <- sum((data - pred)^2)
        ss_tot <- sum((data - mean(data))^2)
        R2_final <- 1 - ss_res / ss_tot
        RMSE_final <- sqrt(mean((data - pred)^2))
      }
    }
  }
  if (!is.finite(tau_final)) {
    if (!is.finite(tau_init) || tau_init <= 0 || tau_init > tau_max) tau_init <- NA_real_
    tau_final <- tau_init
    if (is.finite(tau_final)) {
      pred <- exp(predict(lmfit))
      ss_res <- sum((data - pred)^2)
      ss_tot <- sum((data - mean(data))^2)
      R2_final <- 1 - ss_res / ss_tot
      RMSE_final <- sqrt(mean((data - pred)^2))
    }
  }
  if (return_quality) {
    return(list(tau = tau_final, R2 = R2_final, RMSE = RMSE_final))
  } else {
    return(tau_final)
  }
}


corr_calc <- function(x, y) {
  x <- as.data.frame(x)
  y <- as.numeric(y)
  coref <- numeric(ncol(x))
  for (i in seq_len(ncol(x))) {
    coref[i] <- stats::cor(x[, i], y, use = "pairwise.complete.obs")
  }
  return(coref)
}

get_most_frequent_row <- function(df) {
  combo <- paste(df[[1]], df[[2]], sep = "_")
  tab <- table(combo)
  if (all(tab == 1)) {
    result <- df[which.min(df[[3]]), ]
  } else {
    most_common_combo <- names(tab)[which.max(tab)]
    vals <- strsplit(most_common_combo, "_")[[1]]
    result <- df[df[[1]] == vals[1] & df[[2]] == vals[2], ]
  }
  return(result)
}

calc_acf_rmse <- function(ref, test, frac_acf = 2.5) {
	ok <- is.finite(ref) & is.finite(test)
	ref <- ref[ok]
	test <- test[ok]
	n <- min(length(ref), length(test))
	if (n < 4) return(NA_real_)
	lag.max <- max(floor(n / frac_acf), 5)
	acf_ref  <- acf(ref,  plot = FALSE, lag.max = lag.max)$acf[-1]
	acf_test <- acf(test, plot = FALSE, lag.max = lag.max)$acf[-1]
	m <- min(length(acf_ref), length(acf_test))
	sqrt(mean((acf_ref[1:m] - acf_test[1:m])^2))
}

residual_accuracy_tau <- function(ref, test){
  mat <- cbind(ref, test) %>%
    apply(2, autocorrelation_calc, floor(length(ref) / 2.5))
  tau_results <- apply(mat, 2, function(ac) {
    result <- estimate_tau_lsq(ac, return_quality = TRUE)
    c(tau = result$tau, R2 = result$R2, RMSE = result$RMSE)
  })
  tau_df <- as.data.frame(t(tau_results))
  colnames(tau_df) <- c("tau", "R2", "RMSE")
  tau_df$tau[tau_df$R2 < 0.4 | tau_df$RMSE > 0.5 | !is.finite(tau_df$tau)] <- NA_real_
  tau_final <- tau_df$tau|>as.numeric()
  return(tau_final[1]-tau_final[2])
}

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

evaluate_synthetic_all_metrics <- function(
  Ntim, NpY, id, varphi,L1,L2,t.window, s.window,
  path = NULL
) {
  y <- generate_synthetic_signal(Ntim, NpY, id, varphi)
  dssa  <- decomp_ssa(y$X, frequency = NpY, L1 = L1, L2 = L2)
  dstl  <- decomp_stl(y$X, frequency = NpY, t.window = t.window, s.window = s.window)
  dtrad <- decomp_trad(y$X, frequency = NpY)
  out <- data.frame(
    id = id,
    SSA_trend_cor   = suppressWarnings(cor(dssa$trend,   y$trend,         use = "pairwise.complete.obs")),
    STL_trend_cor   = suppressWarnings(cor(dstl$trend,   y$trend,         use = "pairwise.complete.obs")),
    TRAD_trend_cor  = suppressWarnings(cor(dtrad$trend,  y$trend,         use = "pairwise.complete.obs")),
    SSA_trend_rmse  = calc_rmse(y$trend, dssa$trend),
    STL_trend_rmse  = calc_rmse(y$trend, dstl$trend),
    TRAD_trend_rmse = calc_rmse(y$trend, dtrad$trend),
    SSA_seas_cor    = suppressWarnings(cor(dssa$season,  y$sin_component, use = "pairwise.complete.obs")),
    STL_seas_cor    = suppressWarnings(cor(dstl$season,  y$sin_component, use = "pairwise.complete.obs")),
    TRAD_seas_cor   = suppressWarnings(cor(dtrad$season, y$sin_component, use = "pairwise.complete.obs")),
    SSA_seas_rmse   = calc_rmse(y$sin_component, dssa$season),
    STL_seas_rmse   = calc_rmse(y$sin_component, dstl$season),
    TRAD_seas_rmse  = calc_rmse(y$sin_component, dtrad$season),
    SSA_resid_cor   = suppressWarnings(cor(dssa$residual,  y$residual, use = "pairwise.complete.obs")),
    STL_resid_cor   = suppressWarnings(cor(dstl$residual,  y$residual, use = "pairwise.complete.obs")),
    TRAD_resid_cor  = suppressWarnings(cor(dtrad$residual, y$residual, use = "pairwise.complete.obs")),
    SSA_resid_rmse  = calc_rmse(y$residual, dssa$residual),
    STL_resid_rmse  = calc_rmse(y$residual, dstl$residual),
    TRAD_resid_rmse = calc_rmse(y$residual, dtrad$residual),
    SSA_acf_rmse    = calc_acf_rmse(y$residual, dssa$residual),
    STL_acf_rmse    = calc_acf_rmse(y$residual, dstl$residual),
    TRAD_acf_rmse   = calc_acf_rmse(y$residual, dtrad$residual),
    SSA_delta_tau   = residual_accuracy_tau(dssa$residual,y$residual ),
    STL_delta_tau   = residual_accuracy_tau(dstl$residual,y$residual),
    TRAD_delta_tau  = residual_accuracy_tau(dtrad$residual,y$residual )
  )
  if (!is.null(path)) {
    utils::write.csv(out, path, row.names = FALSE)
  }
  out
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

# =====================real world data=============================

dataproerror <- function(ts,multiplier=1.5){
  Q1 <- quantile(ts, 0.25, na.rm = TRUE)
  Q3 <- quantile(ts, 0.75, na.rm = TRUE)
  IQR_val <- Q3 - Q1
  ts[ts < (Q1 - multiplier * IQR_val) | ts > (Q3 + multiplier * IQR_val)] <- NA
  return(ts)
}

extract_pixel_time_series <- function(rast_data, lon, lat) {
  point <- data.frame(x = lon, y = lat)
  values <- terra::extract(rast_data, point)
  time_series <- as.numeric(values[1, -1])
  if (!is.null(names(rast_data))) {
    names(time_series) <- names(rast_data)
  }
  return(time_series)
}


acfff <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (length(x) < 3) return(NA_real_)
  as.numeric(acf(x, lag.max = 1, plot = FALSE)$acf[2])
}



evaluate_sys_ACF1 <- function(ts, path = NULL) {
  ts <- dataproerror(ts, 1.5)
  if (length(stats::na.omit(ts)) < (length(ts) * 0.75))
    return(rep(NA_real_, 3))
  y <- as.numeric(fill_gap(ts, method = "interp"))
  dssa <- decomp_ssa(y, frequency=NpY,368,69)
  dstl <- decomp_stl(y, frequency=NpY, 367,47)
  dtrad <- decomp_trad(y, frequency=46)
  mat <- cbind(dssa$residual, dstl$residual, dtrad$residual) |>
    apply(2, acfff) |> as.numeric()
  return(mat)
}



evaluate_systest <- function(
  ts,
  frequency = 46,
  ssa_L1 = 368,
  ssa_L2 = 69,
  stl_t.window = 367,
  stl_s.window = 47,
  fill_method = "interp",
  multiplier = 1.5,
  tau_r2_min = 0.4,
  tau_rmse_max = 0.5,
  return_tau_quality = TRUE
) {
  ts <- dataproerror(ts, multiplier = multiplier)
  if (length(stats::na.omit(ts)) < (length(ts) * 0.75)) {
    return(list(
      timeseries = rep(NA_real_, length(ts)),
      residual   = data.frame(SSA = NA_real_, STL = NA_real_, Trad = NA_real_),
      season     = data.frame(SSA = NA_real_, STL = NA_real_, Trad = NA_real_),
      trend      = data.frame(SSA = NA_real_, STL = NA_real_, Trad = NA_real_),
      tau        = c(SSA = NA_real_, STL = NA_real_, Trad = NA_real_),
      tau_quality = data.frame(
        method = c("SSA", "STL", "Trad"),
        tau = NA_real_,
        R2 = NA_real_,
        RMSE = NA_real_
      )
    ))
  }
  y <- as.numeric(fill_gap(ts, method = fill_method))
  Ntim <- length(y)
  dssa  <- decomp_ssa(y, frequency = frequency, L1 = ssa_L1, L2 = ssa_L2)
  dstl  <- decomp_stl(y, frequency = frequency, t.window = stl_t.window, s.window = stl_s.window)
  dtrad <- decomp_trad(y, frequency = frequency)
  res_mat <- cbind(
    SSA  = dssa$residual,
    STL  = dstl$residual,
    Trad = dtrad$residual
  )
  acf_mat <- apply(res_mat, 2, autocorrelation_calc, floor(Ntim / 2.5))
  tau_results <- apply(acf_mat, 2, function(ac) {
    result <- estimate_tau_lsq(ac, return_quality = TRUE)
    c(tau = result$tau, R2 = result$R2, RMSE = result$RMSE)
  })
  tau_df <- as.data.frame(t(tau_results))
  tau_df$method <- rownames(tau_df)
  rownames(tau_df) <- NULL
  tau_df <- tau_df[, c("method", "tau", "R2", "RMSE")]
  tau_final <- tau_df$tau
  bad_tau <- tau_df$R2 < tau_r2_min |
             tau_df$RMSE > tau_rmse_max |
             !is.finite(tau_df$tau)

  tau_final[bad_tau] <- NA_real_
  names(tau_final) <- tau_df$method
  output <- list(
    timeseries = y,
    residual = data.frame(
      SSA  = dssa$residual,
      STL  = dstl$residual,
      Trad = dtrad$residual
    ),
    season = data.frame(
      SSA  = dssa$season,
      STL  = dstl$season,
      Trad = dtrad$season
    ),
    trend = data.frame(
      SSA  = dssa$trend,
      STL  = dstl$trend,
      Trad = dtrad$trend
    ),
    tau = tau_final
  )
  if (return_tau_quality) {
    output$tau_quality <- tau_df
  }
  return(output)
}

evaluate_sys_metrics <- function(ts) {
    ts <- dataproerror(ts, 1.5)
    if (length(na.omit(ts)) < length(ts) * 0.75) return(rep(NA_real_, 8))
    y <- as.numeric(fill_gap(ts,'interp'))
    dssa <- decomp_ssa(y, 46,368,69)$residual
    dstl <- decomp_stl(y, 46, 367,47)$residual
    dtrad <- decomp_trad(y, 46)$residual
    res_mat <- cbind(dssa, dstl, dtrad)
	mat <- res_mat|> apply(2, autocorrelation_calc, floor(length(ts) / 2.5))
	  tau_results <- apply(mat, 2, function(ac) {
		result <- estimate_tau_lsq(ac, return_quality = TRUE)
		c(tau = result$tau, R2 = result$R2, RMSE = result$RMSE)
	  })
	  tau_df <- as.data.frame(t(tau_results))
	  colnames(tau_df) <- c("tau", "R2", "RMSE")
	  tau_df$tau[tau_df$R2 < 0.4 | tau_df$RMSE > 0.5 | !is.finite(tau_df$tau)] <- NA_real_
	  tau_final <- tau_df$tau
    mat1 <- apply(res_mat, 2, acfff)
	out <- c(tau_final,mat1)
    if (length(out) != 8) out <- rep(NA_real_, 8)
    return(out)
}
