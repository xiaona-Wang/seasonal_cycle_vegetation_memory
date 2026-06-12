robin_proj <- "+proj=robin +datum=WGS84"
ssa_rast <- rast("/Volumes/ssa_filter.tif")
q <- values(ssa_rast)[, 1]
q_high <- quantile(q, 0.995, na.rm = TRUE)
q[q > q_high] <- NA
values(ssa_rast) <- q
ssa_df_hist <- data.frame(xyFromCell(ssa_rast, 1:ncell(ssa_rast)), ssa = values(ssa_rast, na.rm = FALSE))

ssa_df_hist <- ssa_df_hist %>%
  mutate(weight = ifelse(is.finite(tau_ssa), cos(y * pi / 180), NA_real_))

wtd_median <- function(x, w){
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (!any(ok)) return(NA_real_)
  x <- x[ok]; w <- w[ok]
  o <- order(x); x <- x[o]; w <- w[o]
  cw <- cumsum(w) / sum(w)
  x[which(cw >= 0.5)[1]]
}
ssa_median <- wtd_median(ssa_df_hist$tau_ssa, ssa_df_hist$weight)
ok <- is.finite(ssa_df_hist$tau_ssa) & is.finite(ssa_df_hist$weight) & ssa_df_hist$weight > 0
den <- density(ssa_df_hist$tau_ssa[ok],
               weights = ssa_df_hist$weight[ok] / sum(ssa_df_hist$weight[ok]),
               from = 0, to = 20, n = 200)
median_density_y <- if (is.finite(ssa_median)) {
  den$y[which.min(abs(den$x - ssa_median))]
} else { NA_real_ }
hist_plot <- ggplot(ssa_df_hist, aes(x = tau_ssa)) +
  geom_histogram(aes(y = ..density.., fill = ..x..),
                 binwidth = 1, color = "grey30", na.rm = TRUE, linewidth = 0.3) +
  scale_fill_stepsn(
    colors = colorRampPalette(brewer.pal(9, "YlGnBu"))(20),
    limits = c(0, 20), breaks = seq(0, 20, by = 1),
    labels = ifelse(seq(0, 20, by = 1) %in% c(0, 5, 10, 15, 20),
        as.character(seq(0, 20, by = 1)), ""),
    oob = scales::squish
  ) +
  geom_vline(xintercept = ssa_median, color = "red", linewidth = 0.4) +
  annotate("text", x = ssa_median, y = median_density_y,
           label = sprintf("%.1f", ssa_median),
           vjust = 0.5, hjust = -0.5, size =3.5, color = "black") +
  scale_x_continuous(limits = c(0, 20), breaks = seq(0, 20, by = 5)) +
  scale_y_continuous(limits = c(0, 0.3), breaks = seq(0, 0.3, by = 0.1)) +
  labs(y = "PDF", x = NULL) +
  theme_minimal(base_size = 8, base_family = "Arial") +
  theme(
    text = element_text(family = "sans"),
    axis.title.y = element_text(size = 10, color = "black"),
    axis.text    = element_text(size = 9, color = "black"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    panel.grid   = element_blank(),
    plot.margin  = margin(2, 2, 2, 2),
    plot.background = element_rect(fill = "white", color = NA),
    legend.position = "none"
  )
  ssa_proj <- project(ssa_rast, robin_proj)
  land <- vect("/Volumes/World_Continents.shp")
  land <- project(land, crs(ssa_proj))
  land_raster <- rasterize(land, ssa_proj, field = 1)
  ssa_df <- as.data.frame(c(ssa_proj, land_raster), xy = TRUE, na.rm = FALSE)
  colnames(ssa_df)[3:4] <- c("ssa", "land_mask")
  ssa_df <- ssa_df %>%
    mutate(category = case_when(
      !is.na(ssa) & ssa == 0 ~ "zero",
      is.na(ssa) & !is.na(land_mask) ~ "land_na",
      !is.na(ssa) ~ "value",
      TRUE ~ "ocean"
    ))
world <- ne_countries(scale = "medium", returnclass = "sf")
world_robin <- st_transform(world, crs = robin_proj)
make_proj_line <- function(lat = NULL, lon = NULL, crs_proj, n = 500) {
  if (!is.null(lat)) {
    lons <- seq(-180, 180, length.out = n)
    lats <- rep(lat, n)
  } else if (!is.null(lon)) {
    lats <- seq(-90, 90, length.out = n)
    lons <- rep(lon, n)
  } else {
    stop("Need to provide either lat or lon.")
  }
  coords <- cbind(lons, lats)
  line <- st_linestring(coords)
  sf_line <- st_sfc(line, crs = 4326)
  st_transform(sf_line, crs = crs_proj)[[1]]
}
grid_lines_robin <- list(
  make_proj_line(lat = 90, crs_proj = robin_proj),
  make_proj_line(lat = -90, crs_proj = robin_proj),
  make_proj_line(lon = 180, crs_proj = robin_proj),
  make_proj_line(lon = -180, crs_proj = robin_proj)
) %>%
  st_sfc(crs = robin_proj)
p_map <- ggplot() +
  geom_raster(data = filter(ssa_df, category == "value"),
    aes(x = x, y = y, fill = ssa)) +
  scale_fill_stepsn(
    colors = colorRampPalette(brewer.pal(9, "YlGnBu"))(20),
    limits = c(0, 20), breaks = seq(0, 20, by = 1),
    labels = ifelse(seq(0, 20, by = 1) %in% c(0, 5, 10, 15, 20),
        as.character(seq(0, 20, by = 1)), ""),
    oob = scales::squish,
    name = expression(tau ~ "(day)")) +
  geom_raster(data = filter(ssa_df, category == "zero"), aes(x = x, y = y), fill = "white") +
  geom_raster(data = filter(ssa_df, category == "land_na"), aes(x = x, y = y), fill = "grey90") +
  geom_sf(data = grid_lines_robin, color = "grey10", linewidth = 0.3, linetype = "solid") +
  geom_sf(data = world_robin, fill = NA, color = "grey15", linewidth = 0.2) +
  coord_sf(crs = robin_proj, expand = FALSE) +
  theme_minimal(base_family = "Arial") +
  theme(
    text = element_text(family = "sans"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b = 2)),
    axis.title = element_blank(),
    axis.text  = element_blank(),
    plot.background = element_rect(fill = "white", color = NA),
    legend.text  = element_text(size = 11, color = "black"),
    legend.title = element_text(size = 13, color = "black"),
    panel.grid   = element_blank(),
    panel.border = element_blank(),
    plot.margin = margin(1, 0, 0, 0),
    legend.position = "right",
    legend.key.height = unit(1.2, "cm"),
    legend.key.width  = unit(0.8, "cm"),
    legend.box.margin = margin(t = -5, r = 0, b = 0, l = -5, unit = "pt")
  ) +
  labs(title = expression("Spatial map of " ~ tau ~ " (SSA)"), x = "", y = "")
p_a <- ggdraw() +
  draw_plot(p_map, x = 0.03, y = 0, width = 1.0, height = 1.0) + 
  draw_plot(hist_plot, x = 0.05, y = 0.1, width = 0.23, height = 0.32)   
 ggsave("/Volumes/F4a_ssa.png",
                    p_a, width = 8, height = 3.6, dpi = 300, bg = "white")

       # Subplot b
       ssa_rast <- rast("/Volumes/ssa_filter.tif")
       ssa_xy <- data.frame(xyFromCell(ssa_rast, 1:ncell(ssa_rast)), ssa = values(ssa_rast, na.rm = FALSE))[,1:2]

       climate <- read.csv("/Volumes/dat_climate19.csv", header = TRUE)
       tcf  <- as.data.frame(read.csv("/Volumes/TCF.csv", header = TRUE))
       sm <- as.data.frame(read.csv("/Volumes/dat_sms.csv", header = TRUE))
       tau <- read.csv("/Volumes/dat_tau.csv", header = TRUE)
       colnames(climate) <- paste0("V", 1:19)
       colnames(sm) <- paste0("sm", 1:3)
       colnames(tcf) <- "tcf"

       num_bins <- 7
       min_x <- min(df_selected$x1);
       max_x <- max(df_selected$x1)
       neg_breaks <- seq(min_x, 0, length.out = 4)
       pos_breaks <- seq(0, max_x, length.out = 4)[-1]
       x1_breaks <- sort(unique(c(neg_breaks, pos_breaks)))
       x1_labels <- paste0(sprintf("%.0f", head(x1_breaks, -1)), "–", sprintf("%.0f", tail(x1_breaks, -1)))

       df <- cbind(climate, tcf, sm, tau, ssa_xy) %>% na.omit() %>% as.data.frame()
       df_selected <- df %>% select(x1 = V1, x2 = V4, y = tau_ssa, lat = y) # MAT vs MAP
       ## mat + map
       x1_breaks <- c(-20, 0, 6, 12, 18, 24, 32)
       x1_labels  <- c("<0","0–6","6–12","12–18","18–24",">24")
       x2_breaks <- c(0, 200, 500, 1000, 2000, 4000, 8500)
       x2_labels  <- c("0–0.2","0.2–0.5","0.5–1.0", "1.0–2.0","2.0–4.0",">4.0")

       df <- cbind(climate, tcf, sm, tau, ssa_xy) %>% na.omit() %>% as.data.frame()
       df_selected <- df %>% select(x1 = V14, x2 = V7, y = tau_ssa, lat = y) # TMP SD vs PRE CV
       # tmp_cv
       x1_breaks <- c(0, 100, 200, 400, 800, 1600, 2400)
       x1_labels <- c("0–1", "1–2", "2–4", "4–8", "8–16", ">16")
       ## pre_cv
       x2_breaks <- c(0, 20, 40, 60, 80, 100, 211)
       x2_labels <- c("0–20", "20–40", "40–60", "60–80", "80–100", ">100")

       df <- cbind(climate, tcf, sm, tau, ssa_xy) %>% na.omit() %>% as.data.frame()
       df_selected <- df %>% select(x1 = sm1, x2 = tcf, y = tau_ssa, lat = y)
       # sm + TCF
       x1_breaks <- seq(min(df_selected$x1), max(df_selected$x1), length.out = num_bins)
       x1_labels <- paste0(sprintf("%.1f", head(x1_breaks, -1)), "–", sprintf("%.1f", tail(x1_breaks, -1)))
       x2_breaks <- seq(min(df_selected$x2), max(df_selected$x2), length.out = num_bins)
       x2_labels <- paste0(sprintf("%.1f", head(x2_breaks, -1)), "–", sprintf("%.1f", tail(x2_breaks, -1)))

       df_selected <- df_selected %>%
         mutate(x1_bin = cut(x1, breaks = x1_breaks, labels = x1_labels, include.lowest = TRUE),
                x2_bin = cut(x2, breaks = x2_breaks, labels = x2_labels, include.lowest = TRUE))
       df_selected$x1_bin <- factor(df_selected$x1_bin, levels = x1_labels)
       df_selected$x2_bin <- factor(df_selected$x2_bin, levels = x2_labels)
        heatmap_data <- df_selected %>%
           mutate(weight = cos(lat * pi / 180)) %>%
           group_by(x1_bin, x2_bin) %>%
           summarise(y_mean = weighted.mean(y, w = weight, na.rm = TRUE),
             count = n(),
             .groups = "drop") %>%
             filter(count >= 20)
       full_grid <- expand_grid(x1_bin = factor(x1_labels, levels = x1_labels),
                                x2_bin = factor(x2_labels, levels = x2_labels))
       heatmap_full <- full_grid %>%
         left_join(heatmap_data, by = c("x1_bin", "x2_bin")) %>%
         mutate(has_data = !is.na(y_mean))
       p_d <- ggplot() +
         geom_tile(data = heatmap_full %>% filter(!has_data), aes(x = x1_bin, y = x2_bin),
                   fill = "grey90", color = "white") +
         geom_tile(data = heatmap_full %>% filter(has_data), aes(x = x1_bin, y = x2_bin, fill = y_mean),
                   color = "white") +
         scale_fill_stepsn(colors = colorRampPalette(brewer.pal(9, "YlGnBu"))(20),
                              limits = c(0,20), breaks = seq(0, 20, by = 1),
                              oob = scales::squish, na.value = "grey90") +

         scale_x_discrete(limits = x1_labels) +
         scale_y_discrete(limits = x2_labels) +
         coord_fixed() +
         #labs(x = "MAT (°C)", y = "MAP (×10³ mm)") +
         #labs(x = "TMP seasonality (°C)", y = "PR seasonality") +
         labs(x = "Mean surface SM (m³/m³)", y = "Tree Cover Fraction") +
         theme_minimal() +
         theme(axis.text.x = element_text(angle = 45, hjust = 1),
               axis.text = element_text(size = 12, color = "black"),
               axis.title = element_text(size = 15, color = "black"),
               panel.grid = element_blank(),
               panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
               panel.background = element_rect(fill = "white", color = NA),
               plot.background = element_rect(fill = "white", color = NA),
               # p_b
               #axis.title.x = element_text(margin = margin(t = -2)),
               #axis.title.y = element_text(margin = margin(r = 5)),
               # p_c
               #axis.title.x = element_text(margin = margin(t = 5)),
               #axis.title.y = element_text(margin = margin(r = 3)),
               # p_d
               axis.title.x = element_text(margin = margin(t = 2)),
               axis.title.y = element_text(margin = margin(r = 3)),
               #plot.margin  = margin(0, 0, 0, 0),
               legend.position = "none")
       ggsave("/Volumes/F4_d.png",
              plot = p_d, width = 4, height = 4,
              dpi = 300)