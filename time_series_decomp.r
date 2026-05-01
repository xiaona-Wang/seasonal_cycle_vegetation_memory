## Time-Series Decomposition

# mean seasonal cycle
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