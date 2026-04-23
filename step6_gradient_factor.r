source("step0_function.r")
##Koppen-Geiger
library(raster)
library(rworldxtra)
library(dplyr)
library(terra)
library(stringr)
period <- '1986-2010'
input_path <- '/net/home/xwang/data/KG/Map_KG-Global/'
output_path <- '/net/home/xwang/decomp/result_new/step6/'
path6 <- "/net/home/xwang/decomp/result_new/step6/"
r <- raster(paste0(input_path, 'KG_', period, '.grd'))
r0 <- r[1:32]; r[1:32] <- seq(1, 32, 1)
r <- ratify(r)
rat <- levels(r)[[1]]
rat$climate <- c('Af', 'Am', 'As', 'Aw', 'BSh', 'BSk', 'BWh', 'BWk',
                 'Cfa', 'Cfb', 'Cfc', 'Csa', 'Csb', 'Csc', 'Cwa', 'Cwb',
                 'Cwc', 'Dfa', 'Dfb', 'Dfc', 'Dfd', 'Dsa', 'Dsb', 'Dsc',
                 'Dsd', 'Dwa', 'Dwb', 'Dwc', 'Dwd', 'EF', 'ET', 'Ocean')
r[1:32] <- r0
levels(r) <- rat
r_template <- raster(extent(r))
res(r_template) <- 0.05
crs(r_template) <- crs(r)
r_resampled <- resample(r, r_template, method = "ngb")
output_tif <- paste0(output_path, "KG_", period, "_resampled_005.tif")
dir.create(output_path, showWarnings = FALSE, recursive = TRUE)
writeRaster(r_resampled, filename = output_tif, format = "GTiff", overwrite = TRUE)
pts <- rasterToPoints(r_resampled, spatial = TRUE)
pts_df <- as.data.frame(pts)
colnames(pts_df)[1:3] <- c("lon", "lat", "KG")
pts_df <- subset(pts_df, KG != 32)
output_csv <- paste0(output_path, "KG_", period, "_resampled_005.csv")
write.csv(pts_df[, c("KG","lon","lat")], file = output_csv, row.names = FALSE)
c <- rast(output_tif)
write.csv(values(c), file.path(output_path, "KB.csv"), row.names = FALSE)

#######PFT
nc_file <- "/net/home/xwang/data/pft/ESACCI-LC-L4-PFT-Map-300m-P1Y-2001-v2.0.8.nc"
pft_stack <- rast(nc_file)
crs(pft_stack) <- "EPSG:4326"
tree_layers <- c("TREES-BE", "TREES-BD", "TREES-NE", "TREES-ND")
tree_stack  <- subset(pft_stack, tree_layers)
tree_cover  <- if (nlyr(tree_stack) > 1) sum(tree_stack, na.rm = TRUE) else tree_stack
tree_cover_fraction <- tree_cover / 100
other_veg_layers <- c("SHRUBS-BE", "SHRUBS-BD", "SHRUBS-NE", "SHRUBS-ND", "GRASS-NAT", "GRASS-MAN")
other_veg_stack  <- subset(pft_stack, other_veg_layers)
other_veg_cover  <- if (nlyr(other_veg_stack) > 1) sum(other_veg_stack, na.rm = TRUE) else other_veg_stack
vegetation_cover_fraction <- (tree_cover + other_veg_cover) / 100
target_res <- 0.05
target_extent <- ext(pft_stack)
target_crs <- crs(pft_stack)
target_template <- rast(ext = target_extent, resolution = target_res, crs = target_crs)
tree_cover_resampled <- resample(tree_cover_fraction, vegetation_cover_fraction, method = "bilinear")
tree_cover_resampled <- resample(tree_cover_resampled, target_template, method = "bilinear")
vegetation_cover_resampled <- resample(vegetation_cover_fraction, target_template, method = "bilinear")
outdir5 <- "/net/home/xwang/decomp/result_new/step6/"
dir.create(outdir5, showWarnings = FALSE, recursive = TRUE)
writeRaster(tree_cover_resampled,       file.path(outdir5, "tree_cover_fraction_2001_0.05deg.tif"), overwrite = TRUE)
writeRaster(vegetation_cover_resampled, file.path(outdir5, "vegetation_cover_fraction_2001_0.05deg.tif"), overwrite = TRUE)
a <- rast(file.path(outdir5, "vegetation_cover_fraction_2001_0.05deg.tif"))
write.csv(values(a), file.path(outdir5, "VCF.csv"), row.names = FALSE)
b <- rast(file.path(outdir5, "tree_cover_fraction_2001_0.05deg.tif"))
write.csv(values(b), file.path(outdir5, "TCF.csv"), row.names = FALSE)

##meteorological data 
path <- '/net/home/xwang/data/worldclim/bio/'
files0 <- list.files(path, full.names = TRUE)
ord <- order(as.numeric(str_extract(basename(files0), "\\d+")))
files <- files0[ord]
lulc <- rast("/net/home/xwang/decomp/result_new/step0/lulc_all.tif")
bi <- lapply(files, function(f) {
  r <- rast(f)
  if (!same.crs(r, lulc)) r <- project(r, lulc, method = "bilinear")
  resample(r, lulc, method = "bilinear")
})
ncell_all <- ncell(lulc)
dat_climate19 <- array(NA_real_, dim = c(ncell_all, length(bi)))
for (i in seq_along(bi)) {
  message("bio var: ", i)
  dat_climate19[, i] <- values(bi[[i]])
}
dir.create(path6, showWarnings = FALSE, recursive = TRUE)
write.csv(dat_climate19, file.path(path6, 'dat_climate19.csv'), row.names = FALSE)

##surface sm
path <- '/net/data/GLEAM/data/v4.1a/monthly/'
files <- list.files(path)
files <- files[12]
lulc <- rast("/net/home/xwang/decomp/result_new/step0/lulc_all.tif")
dat_sm <- NULL
for (i in seq_along(files)) {
  message("GLEAM year dir: ", files[i])
  y <- list.files(file.path(path, files[i]), full.names = TRUE)[-c(1:21)]
  ncell_all <- ncell(lulc)
  test_r <- rast(y[1]) %>% project(lulc, method = "bilinear")
  L <- nlyr(test_r)
  datall <- array(NA_real_, dim = c(ncell_all, length(y) * L))
  for (j in seq_along(y)) {
    message("  file ", j, "/", length(y))
    r <- rast(y[j]) %>% project(lulc, method = "bilinear")
    dats <- array(NA_real_, dim = c(ncell_all, nlyr(r)))
    for (k in 1:nlyr(r)) {
      dats[, k] <- values(r[[k]])
    }
    idx <- ((j - 1) * nlyr(r) + 1):(j * nlyr(r))
    datall[, idx] <- dats
  }
  funcv <- function(x) sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE)
  datsm_result1 <- apply(datall, 1, mean, na.rm = TRUE)
  nyear <- ncol(datall) / 12
  dat_yearmean <- array(NA_real_, dim = c(ncell_all, nyear))
  for (ii in 1:nyear) {
    idx <- ((ii - 1) * 12 + 1):(ii * 12)
    dat_yearmean[, ii] <- apply(datall[, idx, drop = FALSE], 1, mean, na.rm = TRUE)
  }
  datsm_result2 <- apply(dat_yearmean, 1, funcv)
  month_means <- sapply(1:12, function(m) {
    idx <- seq(m, ncol(datall), by = 12)
    apply(datall[, idx, drop = FALSE], 1, mean, na.rm = TRUE)
  })
  datsm_result3 <- apply(month_means, 1, funcv)
  dat_sm <- cbind(datsm_result1, datsm_result2, datsm_result3)
}
write.csv(dat_sm, file.path(path6, 'dat_sms.csv'), row.names = FALSE)
