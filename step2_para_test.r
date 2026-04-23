source("step0_function.r")
library(dplyr)
library(doParallel)
library(parallel)
library(foreach)

path2 <- '/Volumes/STUDY1/topic1/output/step2/'
cl <- makeCluster(10)
registerDoParallel(cl)
datssa_layered <- foreach(
  x = 1:100,
  .combine = "rbind",
  .packages = c("zoo", "dplyr", "doParallel", "Rssa", "forecast")
) %dopar% {
  ans <- find_best_ssa_params(
    Ntim = 46 * 23,
    NpY = 46,
    id = x,
    varphi = 0.6
  )
  if (is.null(ans)) {
    data.frame(method = "SSA", L1 = NA_real_, L2 = NA_real_)
  } else {
    ans$best[, c("method", "L1", "L2")]
  }
}
stopCluster(cl)
write.csv(datssa_layered, paste0(path2, 'datssa_layered.csv'), row.names = FALSE)

################
cl <- makeCluster(10)
registerDoParallel(cl)
datstl_layered <- foreach(
  x = 1:100,
  .combine = "rbind",
  .packages = c("zoo", "dplyr", "doParallel", "Rssa", "forecast")
) %dopar% {
  ans <- find_best_stl_params(
    Ntim = 46 * 23,
    NpY = 46,
    id = x,
    varphi = 0.6
  )
  if (is.null(ans)) {
    data.frame(method = "STL", t.window = NA_real_, s.window = NA_real_)
  } else {
    ans$best[, c("method", "t.window", "s.window")]
  }
}
stopCluster(cl)
write.csv(datstl_layered, paste0(path2, 'datstl_layered.csv'), row.names = FALSE)

get_most_frequent_params <- function(value) {
	names(value)[2:3] <- c('trend','season')
  dat <- value |>
    dplyr::filter(is.finite(trend), is.finite(season))
  if (nrow(dat) == 0) {
    return(NULL)
  }
  freq_tab <- dat |>
    dplyr::count(trend, season, name = "freq") |>
    dplyr::mutate(prop = freq / sum(freq)) |>
    dplyr::arrange(dplyr::desc(freq), trend, season)
  list(
    best = freq_tab[1, ],
    all = freq_tab
  )
}
################
ssa_freq <- get_most_frequent_params(datssa_layered)
ssa_freq$best   
ssa_freq$all   
stl_freq <- get_most_frequent_params(datstl_layered)
stl_freq$best
stl_freq$all