source("step0_function.r")
library(parallel)
library(doParallel)
library(foreach)
library(stringr)
path4  <- '/Volumes/STUDY1/topic1/output/step4/'
.export_funcs <- c(
    "calc_rmse", "calc_acf_rmse", "residual_accuracy_tau",
    "evaluate_synthetic_all_metrics",
    "generate_synthetic_signal",
    "decomp_ssa", "decomp_stl", "decomp_trad",
    "fill_gap", "fill_gap_auto",
    "autocorrelation_calc", "estimate_tau_lsq"
)
.packages_common <- c("zoo","dplyr","doParallel","Rssa","tseries","urca","forecast")
### example 1 (TREND)
dir.create(paste0(path4,'Example1'), showWarnings = FALSE)
path41 <- paste0(path4,'Example1/')
cl <- makeCluster(10)
registerDoParallel(cl)
trend_scale_seq <- c(0.14, 0.16, 0.18, 0.20, 0.22)
for (i in 1:100) {
  print(i)
  invisible(foreach(
    x = 1:length(trend_scale_seq),
    .export   = .export_funcs,
    .packages = .packages_common
  ) %dopar% {
    evaluate_synthetic_tau_test(Ntim=46*23, NpY=46, id=i, varphi=0.6, trend_scale = trend_scale_seq[x],
                           path=paste0(path41, formatC(i,width=3,flag='0'),'_', sprintf("%05d_res.csv", x)))
  })
}
stopCluster(cl)
files   <- list.files(path41)
filesub <- substr(files, 5, str_length(files) - 8)
datmedian <- datmean <- datsd <- array(data = NA_real_, dim = c(length(unique(filesub)), 4))
for (i in seq_along(unique(filesub))) {
  print(i)
  y   <- which(filesub == unique(filesub)[i])
  dat <- array(data = NA_real_, dim = c(length(y), 4))
  for (j in seq_along(y)) {
    dat[j, ] <- read.csv(paste0(path41, files[y[j]]), header = TRUE)[,1]
  }
  datmedian[i, ] <- apply(dat, 2, median, na.rm = TRUE)
  datmean[i, ]   <- apply(dat, 2, mean,   na.rm = TRUE)
  datsd[i, ]     <- apply(dat, 2, sd,     na.rm = TRUE)
}
write.csv(datmedian, paste0(path4,'Example1_median.csv'), row.names = FALSE)
write.csv(datmean,   paste0(path4,'Example1_mean.csv'),   row.names = FALSE)
write.csv(datsd,     paste0(path4,'Example1_sd.csv'),     row.names = FALSE)
### example 2 (SEASONAL)
dir.create(paste0(path4,'Example2'), showWarnings = FALSE)
path42 <- paste0(path4,'Example2/')
amp_mod_seq <- c(0.05, 0.1, 0.3, 0.6, 0.8)  
cl <- makeCluster(10)
registerDoParallel(cl)
for (i in 1:100) {
  print(i)
  invisible(foreach(
    x = 1:length(amp_mod_seq),  
    .export   = .export_funcs,
    .packages = .packages_common
  ) %dopar% {
     evaluate_synthetic_tau_test(Ntim=46*23, NpY=46, id=i, varphi=0.6,
								amp_mod = amp_mod_seq[x], 
                           path=paste0(path42, formatC(i,width=3,flag='0'),'_', sprintf("%05d_res.csv", x)))
  })
}
stopCluster(cl)
files   <- list.files(path42)
filesub <- substr(files, 5, str_length(files) - 8)
datmedian <- datmean <- datsd <- array(data = NA_real_, dim = c(length(unique(filesub)), 4))
for (i in seq_along(unique(filesub))) {
  print(i)
  y   <- which(filesub == unique(filesub)[i])
  dat <- array(data = NA_real_, dim = c(length(y), 4))
  for (j in seq_along(y)) {
    dat[j, ] <- read.csv(paste0(path42, files[y[j]]), header = TRUE)[,1]
  }
  datmedian[i, ] <- apply(dat, 2, median, na.rm = TRUE)
  datmean[i, ]   <- apply(dat, 2, mean,   na.rm = TRUE)
  datsd[i, ]     <- apply(dat, 2, sd,     na.rm = TRUE)
}
write.csv(datmedian, paste0(path4,'Example2_median.csv'), row.names = FALSE)
write.csv(datmean,   paste0(path4,'Example2_mean.csv'),   row.names = FALSE)
write.csv(datsd,     paste0(path4,'Example2_sd.csv'),     row.names = FALSE)
### example 3 (VARPHI)
dir.create(paste0(path4,'Example3'), showWarnings = FALSE)
path43 <- paste0(path4,'Example3/')
cl <- makeCluster(10)
registerDoParallel(cl)
varphi <- c(0.4, 0.5, 0.6, 0.7, 0.8)
for (i in 1:100) {
  print(i)
  invisible(foreach(
    x = 1:length(varphi),
    .export   = .export_funcs,
    .packages = .packages_common
  ) %dopar% {
   evaluate_synthetic_tau_test(Ntim=46*23, NpY=46, id=i, varphi= varphi[x],
                           path=paste0(path43, formatC(i,width=3,flag='0'),'_', sprintf("%05d_res.csv", x)))
  })
}
stopCluster(cl)
files   <- list.files(path43)
filesub <- substr(files, 5, str_length(files) - 8)
datmedian <- datmean <- datsd <- array(data = NA_real_, dim = c(length(unique(filesub)), 4))
for (i in seq_along(unique(filesub))) {
  print(i)
  y   <- which(filesub == unique(filesub)[i])
  dat <- array(data = NA_real_, dim = c(length(y), 4))
  for (j in seq_along(y)) {
    dat[j, ] <- read.csv(paste0(path43, files[y[j]]), header = TRUE)[,1]
  }
  datmedian[i, ] <- apply(dat, 2, median, na.rm = TRUE)
  datmean[i, ]   <- apply(dat, 2, mean,   na.rm = TRUE)
  datsd[i, ]     <- apply(dat, 2, sd,     na.rm = TRUE)
}
write.csv(datmedian, paste0(path4,'Example3_median.csv'), row.names = FALSE)
write.csv(datmean,   paste0(path4,'Example3_mean.csv'),   row.names = FALSE)
write.csv(datsd,     paste0(path4,'Example3_sd.csv'),     row.names = FALSE)
### example 4 (RESOLUTION)
dir.create(paste0(path4,'Example4'), showWarnings = FALSE)
path44 <- paste0(path4,'Example4/')
cl <- makeCluster(10)
registerDoParallel(cl)
Temporalresolution <- c(24, 36, 46, 73, 180)
for (i in 1:100) {
  print(i)
  invisible(foreach(
    x = 1:length(Temporalresolution),
    .export   = .export_funcs,
    .packages = .packages_common
  ) %dopar% {
    evaluate_synthetic_tau_test(Ntim=Temporalresolution[x]*23, NpY=Temporalresolution[x],  id=i,   varphi= 0.6,
							L1=floor((368/46)*Temporalresolution[x]), L2=floor((69/46)*Temporalresolution[x]),
							t.window=floor((367/46)*Temporalresolution[x]), s.window=floor((47/46)*Temporalresolution[x]),
                            path = paste0(path44, formatC(i,width=3,flag='0'),'_', sprintf("%05d_res.csv", x)))
  })
}
stopCluster(cl)
files   <- list.files(path44)
filesub <- substr(files, 5, str_length(files) - 8)
datmedian <- datmean <- datsd <- array(data = NA_real_, dim = c(length(unique(filesub)), 4))
for (i in seq_along(unique(filesub))) {
  print(i)
  y   <- which(filesub == unique(filesub)[i])
  dat <- array(data = NA_real_, dim = c(length(y), 4))
  for (j in seq_along(y)) {
    dat[j, ] <- read.csv(paste0(path44, files[y[j]]), header = TRUE)[,1]
  }
  datmedian[i, ] <- apply(dat, 2, median, na.rm = TRUE)
  datmean[i, ]   <- apply(dat, 2, mean,   na.rm = TRUE)
  datsd[i, ]     <- apply(dat, 2, sd,     na.rm = TRUE)
}
write.csv(datmedian, paste0(path4,'Example4_median.csv'), row.names = FALSE)
write.csv(datmean,   paste0(path4,'Example4_mean.csv'),   row.names = FALSE)
write.csv(datsd,     paste0(path4,'Example4_sd.csv'),     row.names = FALSE)
### example 5 (LENGTH)
dir.create(paste0(path4,'Example5'), showWarnings = FALSE)
path45 <- paste0(path4,'Example5/')
cl <- makeCluster(10)
registerDoParallel(cl)
datalength <- seq(10, 50, 10)
for (i in 1:100) {
  print(i)
  invisible(foreach(
    x = 1:length(datalength),
    .export   = .export_funcs,
    .packages = .packages_common
  ) %dopar% {
    evaluate_synthetic_tau_test(Ntim=46*datalength[x], NpY=46,  id=i, varphi= 0.6,
                         path=  paste0(path45, formatC(i,width=3,flag='0'),'_', sprintf("%05d_res.csv", x)))
  })
}
stopCluster(cl)
files   <- list.files(path45)
filesub <- substr(files, 5, str_length(files) - 8)
datmedian <- datmean <- datsd <- array(data = NA_real_, dim = c(length(unique(filesub)), 4))
for (i in seq_along(unique(filesub))) {
  print(i)
  y   <- which(filesub == unique(filesub)[i])
  dat <- array(data = NA_real_, dim = c(length(y), 4))
  for (j in seq_along(y)) {
    dat[j, ] <- read.csv(paste0(path45, files[y[j]]), header = TRUE)[,1]
  }
  datmedian[i, ] <- apply(dat, 2, median, na.rm = TRUE)
  datmean[i, ]   <- apply(dat, 2, mean,   na.rm = TRUE)
  datsd[i, ]     <- apply(dat, 2, sd,     na.rm = TRUE)
}
write.csv(datmedian, paste0(path4,'Example5_median.csv'), row.names = FALSE)
write.csv(datmean,   paste0(path4,'Example5_mean.csv'),   row.names = FALSE)
write.csv(datsd,     paste0(path4,'Example5_sd.csv'),     row.names = FALSE)
###### example 6 (NOISE)
dir.create(paste0(path4,'Example6'), showWarnings = FALSE)
path46 <- paste0(path4,'Example6/')
cl <- makeCluster(10)
registerDoParallel(cl)
w_n_sd <- c(0, 0.005, 0.01, 0.02, 0.05)
for (i in 1:100) {
  print(i)
  invisible(foreach(
    x = 1:length(w_n_sd),
    .export   = .export_funcs,
    .packages = .packages_common
  ) %dopar% {
   evaluate_synthetic_tau_test(Ntim=46*23, NpY=46, id=i, varphi=0.6, white_noise_sd = w_n_sd[x],
                          path= paste0(path46, formatC(i,width=3,flag='0'),'_', sprintf("%05d_res.csv", x)))
  })
}
stopCluster(cl)
files   <- list.files(path46)
filesub <- substr(files, 5, str_length(files) - 8)
datmedian <- datmean <- datsd <- array(data = NA_real_, dim = c(length(unique(filesub)), 4))
for (i in seq_along(unique(filesub))) {
  print(i)
  y   <- which(filesub == unique(filesub)[i])
  dat <- array(data = NA_real_, dim = c(length(y), 4))
  for (j in seq_along(y)) {
    dat[j, ] <- read.csv(paste0(path46, files[y[j]]), header = TRUE)[,1]
  }
  datmedian[i, ] <- apply(dat, 2, median, na.rm = TRUE)
  datmean[i, ]   <- apply(dat, 2, mean,   na.rm = TRUE)
  datsd[i, ]     <- apply(dat, 2, sd,     na.rm = TRUE)
}
write.csv(datmedian, paste0(path4,'Example6_median.csv'), row.names = FALSE)
write.csv(datmean,   paste0(path4,'Example6_mean.csv'),   row.names = FALSE)
write.csv(datsd,     paste0(path4,'Example6_sd.csv'),     row.names = FALSE)