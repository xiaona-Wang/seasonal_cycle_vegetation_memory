source("step0_function.r")
library(terra)
library(parallel)
path_in   <- "/net/home/xwang/decomp/result_new/step1/"
path_out  <- "/net/home/xwang/decomp/result_new/step5/"
lulc      <- "/net/home/xwang/decomp/result_new/step0/lulc_all.tif"
continents_shp <- "/net/home/xwang/data/world/World_Continents.shp"
path6 <- "/net/home/xwang/decomp/result_new/step6/"
if (!dir.exists(path_out)) dir.create(path_out, recursive = TRUE)
world_continent <- terra::vect(continents_shp)
world_continent <- terra::aggregate(world_continent)
r_lulc <- terra::rast(lulc)
stacktif <- list.files(path_in, "\\.tif$", full.names = TRUE)[-c(1:46)]
r_stack  <- terra::rast(stacktif)
if (!terra::same.crs(r_lulc, r_stack)) {
 cat("Projecting LULC...\n")
 r_lulc <- terra::project(r_lulc, r_stack, method = "near")
}
if (!terra::same.crs(world_continent, r_stack)) {
 world_continent <- terra::project(world_continent, r_stack)
}
r_mask <- ifel(r_lulc == 1, 1, NA)
r_stack_masked <- mask(r_stack, r_mask)
qts_australia <- extract_pixel_time_series(r_stack_masked, 142.075, -33.025)|>as.numeric()|>evaluate_systest()
write.csv(data.frame(timeseries=qts_australia$timeseries,residual=qts_australia$residual,
            seasonal=qts_australia$season,trend=qts_australia$trend),
            	'/net/home/xwang/decomp/result_new/step5/hotspot/australia.csv',row.names=FALSE)
for(i in 1:3600){
	write.csv(r_stack_masked[i,,],paste0(path_out,formatC(i,width=4,flag='0'),'.csv'),row.names=F)
}
files <- list.files(path6, pattern = "\\.csv$", full.names = TRUE)
process_file <- function(f) {
  df <- read.csv(f,header=T)
  dat <- apply(df, 1, evaluate_sys_metrics)
  out_file <- file.path(
    path_out,
    paste0(basename(tools::file_path_sans_ext(f)), "_out.csv")
  )
  write.csv(dat, out_file, row.names = FALSE)
}
cl <- makeCluster(20)
clusterExport(cl, c("evaluate_sys_metrics", "process_file", "path_out","dataproerror",
					"decomp_ssa","decomp_stl","decomp_trad","acfff",
					"autocorrelation_calc","estimate_tau_lsq",
					"evaluate_synthetic_all_metrics",
					"calc_rmse", "calc_acf_rmse", "residual_accuracy_tau",
                    "fill_gap_auto","fill_gap","mean_cycle_calc"))
clusterEvalQ(cl, {
  library(Rssa)
  library(forecast)
})
res <- parLapply(cl, files, process_file)
stopCluster(cl)