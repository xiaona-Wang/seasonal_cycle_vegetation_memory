## Memory Estimate

# autocorrelation function calculation
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

## tau estimate calculation
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