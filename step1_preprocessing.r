
source("step0_function.r")
### data preprocessing
library(reticulate)
library(dplyr)
library(raster)
library(terra)
use_condaenv("xiaona", required = TRUE)
py_config()
zarr <- import("zarr")
path    <- '/net/data/MODIS/MCD43C4/kndvi-8D-0.05deg-46x600x600/data/'
outpath <- '/net/home/xwang/decomp/result_new/step1/'
files   <- list.files(path, pattern = "\\.zarr$", full.names = TRUE)
dattif <- raster(ncol = 600, nrow = 600)
for (i in seq_along(files)) {
  cat("zarr file:", i, " / ", length(files), "\n")
  dat <- zarr$open(files[i], mode = "r")
  dataset <- dat$get("kNDVI")
  n_t <- dataset$shape[[1]]
  for (j in seq_len(n_t)) {
    arr <- reticulate::py_to_r(dataset[j, , ])
    values(dattif) <- flatten_2d_to_1d(arr)
    writeRaster(
      dattif,
      filename  = file.path(outpath, paste0(1999 + i, formatC(j, width = 2, flag = '0'), '.tif')),
      overwrite = TRUE
    )
  }
}


path_lulc_dir <- "/net/home/xwang/data/MCD12C1_lulc_tif/"
path_ts_dir   <- "/net/home/xwang/decomp/result_new/step1/"
outpath       <- "/net/home/xwang/decomp/result_new/step0/"
dir.create(outpath, showWarnings = FALSE, recursive = TRUE)
cores         <- 1
tempdir_fast  <- "/fast/tmp"
memfrac       <- 0.8
wopt_mask <- list(
  datatype = "INT1U",
  gdal = c("COMPRESS=LZW", "NUM_THREADS=ALL_CPUS", "BIGTIFF=YES")
)
if (!dir.exists(tempdir_fast)) dir.create(tempdir_fast, recursive = TRUE)
terraOptions(tempdir = tempdir_fast, memfrac = memfrac)
lulc_files <- list.files(path_lulc_dir, full.names = TRUE)
ts_files   <- list.files(path_ts_dir, pattern = "\\.tif$", full.names = TRUE)
r_lulc   <- rast(lulc_files)
ts_stack <- rast(ts_files)
if (!compareGeom(r_lulc[[1]], ts_stack[[1]], stopOnError = FALSE)) {
  r_lulc <- project(r_lulc, ts_stack[[1]], method = "near")
}
vege_file <- file.path(outpath, "lulc_all.tif")
vege <- detect_vege_pixels(
  tif            = r_lulc,
  veg_codes      = 1:11,
  min_valid_frac = 0.75,
  cores          = cores,
  out            = vege_file,
  wopt           = wopt_mask,
  .set_options   = TRUE,
  tempdir        = tempdir_fast,
  memfrac        = memfrac
)

ts_stack_vege <- mask(ts_stack, vege, maskvalue = 0, updatevalue = NA)
flat_file <- file.path(outpath, "lulc_flat.tif")
flat_mask <- detect_flat_pixels(
  tif            = ts_stack_vege,
  min_valid_frac = 0.75,
  thr_mean_high  = 0.5,
  thr_mean_low   = 0.03,
  thr_sd         = 0.0001,
  thr_dx         = 0.0001,
  run_len        = 13 - 1,
  cores          = cores,
  out            = flat_file,
  wopt           = wopt_mask,
  .set_options   = TRUE,
  tempdir        = tempdir_fast,
  memfrac        = memfrac
)
analysis_mask <- ifel(vege == 1 & flat_mask == 0, 1, NA)
analysis_mask_file <- file.path(outpath, "analysis_mask.tif")
writeRaster(
  analysis_mask,
  filename = analysis_mask_file,
  overwrite = TRUE,
  wopt = wopt_mask
)
ts_stack_final <- mask(ts_stack, analysis_mask)


path_lulc_dir <- "/net/home/xwang/data/MCD12C1_lulc_tif/"
outpath       <- "/net/home/xwang/decomp/result_new/step0/"

wopt_mask <- list(
  datatype = "INT1U",
  gdal = c("COMPRESS=LZW", "NUM_THREADS=ALL_CPUS", "BIGTIFF=YES")
)
lulc_files <- sort(list.files(path_lulc_dir, pattern = "\\.tif$", full.names = TRUE))
r_lulc <- rast(lulc_files)
analysis_mask <- rast(file.path(outpath, "analysis_mask.tif"))
if (!compareGeom(r_lulc[[1]], analysis_mask, stopOnError = FALSE)) {
  r_lulc <- project(r_lulc, analysis_mask, method = "near")
}
veg_codes <- 1:11
r_lulc_veg <- ifel(r_lulc %in% veg_codes, r_lulc, NA)
lulc_mode <- modal(r_lulc_veg, na.rm = TRUE)
lulc_analysis <- mask(lulc_mode, analysis_mask)
lulc_analysis_file <- file.path(outpath, "analysis_mask_lulc_mode.tif")
writeRaster(
  lulc_analysis,
  filename = lulc_analysis_file,
  overwrite = TRUE,
  wopt = wopt_mask
)