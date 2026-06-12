robin_proj <- "+proj=robin +datum=WGS84"
#ssa_raw <- rast("/Volumes/stl_dif_filter.tif")
#ssa_raw <- rast("/Volumes/trad_dif_filter.tif")
ssa_raw <- rast("/Volumes/trad_stl_filter.tif")
q <- values(ssa_raw)[, 1]
q_low  <- quantile(q, 0.01, na.rm = TRUE)
q_high <- quantile(q, 0.99, na.rm = TRUE)
q[q < q_low | q > q_high] <- NA
values(ssa_raw) <- q
ssa_df_hist <- as.data.frame(ssa_raw, xy = TRUE, na.rm = FALSE)
colnames(ssa_df_hist)[3] <- "ssa"
ssa_df_hist <- ssa_df_hist %>%
  mutate(weight = ifelse(is.finite(ssa), cos(y * pi / 180), NA_real_))
wtd_median <- function(x, w){
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (!any(ok)) return(NA_real_)
  x <- x[ok]; w <- w[ok]
  o <- order(x); x <- x[o]; w <- w[o]
  cw <- cumsum(w) / sum(w)
  x[which(cw >= 0.5)[1]]
}
ssa_median <- wtd_median(ssa_df_hist$ssa, ssa_df_hist$weight)
ok <- is.finite(ssa_df_hist$ssa) & is.finite(ssa_df_hist$weight) & ssa_df_hist$weight > 0
den <- density(ssa_df_hist$ssa[ok],
              weights = ssa_df_hist$weight[ok] / sum(ssa_df_hist$weight[ok]),
              from = -40, to = 40, n = 600)
median_density_y <- if (is.finite(ssa_median)) {
  den$y[which.min(abs(den$x - ssa_median))]
} else { NA_real_ }
y_text <- min(median_density_y, 0.03)
hist_plot <- ggplot(ssa_df_hist, aes(x = ssa, weight = weight)) +
  geom_histogram(
    aes(y = ..density.., fill = ..x..),
    binwidth = 3, color = "grey30", na.rm = TRUE, linewidth = 0.3
  ) +
  scale_fill_stepsn(
    colors = colorRampPalette(rev(brewer.pal(11, "RdBu")))(20),
    limits = c(-40, 40),
    breaks = seq(-40, 40, by = 4),
    oob = scales::squish
  ) +
  geom_vline(xintercept = ssa_median, color = "red", linewidth = 0.4) +
  annotate(
    "text",
    x = ssa_median,
    y = y_text,
    label = sprintf("%.1f", ssa_median),
    ## trad_stl
    vjust = -0.8, hjust = 2.0,
    ## trad_ssa
    #vjust = -0.5, hjust = 2.5,
    ## stl_ssa
    #vjust = 0, hjust = 2.0,
    size = 2.5,
    color = "black"
  ) +
  #scale_y_continuous(
  #  #limits = c(0, 0.06),
  #  breaks = seq(0, 0.06, by = 0.03)
  #) +
  scale_y_continuous(
  breaks = scales::pretty_breaks(n = 3)
)+
  scale_x_continuous(
    limits = c(-40, 40),
    breaks = seq(-40, 40, by = 20)
  ) +
  labs(y = "PDF", x = NULL) +
  theme_minimal(base_size = 7, base_family = "Arial") +
  theme(
    axis.title.y = element_text(size = 6.5, color = "black"),
    axis.text = element_text(size = 6.5, color = "black"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    panel.grid = element_blank(),
    plot.margin = margin(2, 2, 2, 2),
    plot.background = element_rect(fill = "white", color = NA),
    legend.position = "none"
  )
ssa <- project(ssa_raw, robin_proj)
land <- vect("/Volumes/World_Continents.shp")
land <- project(land, crs(ssa))
land_raster <- rasterize(land, ssa, field = 1)
ssa_df <- as.data.frame(c(ssa, land_raster), xy = TRUE, na.rm = FALSE)
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
p <- ggplot() +
  geom_raster(data = filter(ssa_df, category == "value"),
    aes(x = x, y = y, fill = ssa)) +
  scale_fill_stepsn(
    colors = colorRampPalette(rev(brewer.pal(11, "RdBu")))(20),
    limits = c(-40, 40), breaks = seq(-40, 40, by = 4),
    labels = ifelse(seq(-40, 40, by = 4) %in% c(-40, -20, 0, 20, 40),
                    as.character(seq(-40, 40, by = 4)), ""),
    oob = scales::squish,
    name = expression(Delta ~ tau ~ "(day)")
  ) +
  geom_raster(data = filter(ssa_df, category == "zero"),aes(x = x, y = y), fill = "white") +
  geom_raster(data = filter(ssa_df, category == "land_na"),aes(x = x, y = y), fill = "grey90") +
  geom_sf(data = grid_lines_robin, color = "grey10", linewidth = 0.3, linetype = "solid") +
  geom_sf(data = world_robin, fill = NA, color = "grey15", linewidth = 0.2) +

  coord_sf(crs = robin_proj, expand = FALSE) +
  theme_minimal(base_family = "Arial") +
  theme(
    text = element_text(family = "sans"),
    plot.title = element_text(size = 10, face = "bold", hjust = 0.5, margin = margin(b = 2)),
    axis.title = element_blank(),
    axis.text = element_blank(),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    plot.margin = margin(1, 0, 0, 0),
    legend.position = "bottom",
    legend.text = element_text(size = 8, color = "black"),
    legend.title = element_text(size = 9, color = "black"),
    legend.key.height = unit(0.4, "cm"),
    legend.key.width  = unit(1.6, "cm"),
    legend.box.margin = margin(t = -5, unit = "pt")
  ) +
  labs(
    #title = expression("Spatial map of "* Delta ~ tau ~ " (STL - SSA)"),
    #title = expression("Spatial map of "* Delta ~ tau ~ " (Trad - SSA)"),
    title = expression("Spatial map of "* Delta ~ tau ~ " (Trad - STL)"),
    x = "", y = ""
  )
final_plot <- ggdraw(p) +
  draw_plot(hist_plot, x = 0.03, y = 0.25, width = 0.26, height = 0.30)
ggsave("/Volumes/F3_trad_stl.png", final_plot,
       width = 5, height = 3, dpi = 300,
       bg = "white")

# subplot b_lulc
path6 <- "/Volumes/"
lulc_type <- rast("/Volumes/analysis_mask_lulc_mode.tif")
ssa_rast <- rast(paste0(path6, "ssa_filter.tif"))
tau <- read.csv(paste0(path6, "dat_tau.csv"), header = TRUE)
ssa_df <- data.frame(xyFromCell(ssa_rast, 1:ncell(ssa_rast)), ssa = values(ssa_rast, na.rm = FALSE))
ssa_df <- ssa_df %>%
  mutate(weight = ifelse(is.finite(tau_ssa), cos(y * pi / 180), NA_real_))
lulc_df <- as.data.frame(lulc_type, xy = FALSE, na.rm = FALSE)
colnames(lulc_df) <- "lulc_type"
df <- cbind(lulc_df, tau, ssa_df[,c(1,2,4)]) %>%
  na.omit()
df_sele <- df %>%
  mutate(
  x1_group = case_when(
      lulc_type %in% 1:2 ~ "evergreen",
      lulc_type %in% 3:5 ~ "deciduous",
      lulc_type %in% 6:7 ~ "shrublands",
      lulc_type %in% 8:9 ~ "savannas",
      lulc_type == 10 ~ "grasslands",
      lulc_type == 11 ~ "other",
      TRUE ~ NA_character_
  ),
    y1 = tau_trad.1,   
    y2 =  tau_stl.1    
  ) %>%
  filter(!is.na(x1_group), x1_group != "other")
q_y1 <- quantile(df_sele$y1, probs = c(0.01, 0.85), na.rm = TRUE)
q_y2 <- quantile(df_sele$y2, probs = c(0.01, 0.85), na.rm = TRUE)
df_y1 <- df_sele %>%
  filter(y1 >= q_y1[1], y1 <= q_y1[2]) %>%
  transmute(x1_group, y = y1, type = "Trad", weight)
df_y2 <- df_sele %>%
  filter(y2 >= q_y2[1], y2 <= q_y2[2]) %>%
  transmute(x1_group, y = y2, type = "STL", weight)
df_long <- bind_rows(df_y1, df_y2) %>%
  mutate(
    type = factor(type, levels = c("Trad", "STL")),
    x1_group = factor(x1_group, levels = c("evergreen","deciduous","shrublands","savannas", "grasslands"))
  )
df_long$x1_group <- recode(df_long$x1_group,
                           "evergreen" = "EF",
                           "deciduous" = "DF",
                           "shrublands" = "SH",
                           "savannas" = "SA",
                           "grasslands" = "GR")
wtd_median <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (!any(ok)) return(NA_real_)
  x <- x[ok]; w <- w[ok]
  o <- order(x); x <- x[o]; w <- w[o]
  cw <- cumsum(w) / sum(w)
  x[which(cw >= 0.5)[1]]
}
medians <- df_long %>%
  group_by(x1_group, type) %>%
  summarise(med = wtd_median(y, weight), .groups = "drop")
pw_methods <- df_long %>%
  group_by(x1_group) %>%
  pairwise_wilcox_test(y ~ type, p.adjust.method = "BH") %>%
  ungroup()
pw_groups <- df_long %>%
  group_by(type) %>%
  pairwise_wilcox_test(y ~ x1_group, p.adjust.method = "BH") %>%
  ungroup()
 letters_methods <- pw_methods %>%
    group_by(x1_group) %>%
    group_modify(~{
      nm <- paste(.x$group1, .x$group2, sep = "-")
      sig <- setNames(.x$p.adj <= 0.05, nm)
      L <- multcompView::multcompLetters(sig)$Letters
      expected_order <- c("Trad", "STL")
      L <- L[intersect(expected_order, names(L))]
      tibble(type = names(L), letter_upper = toupper(unname(L)))
    }) %>%
    ungroup()
letters_groups <- pw_groups %>%
  group_by(type) %>%
  group_modify(~{
    nm <- paste(.x$group1, .x$group2, sep = "-")
    sig <- setNames(.x$p.adj <= 0.05, nm)
    L <- multcompView::multcompLetters(sig)$Letters
    tibble(x1_group = names(L), letter_lower = tolower(unname(L)))
  }) %>%
  ungroup()
letters_df <- medians %>%
  mutate(letter_upper = case_when(type == "Trad" ~ "A", type == "STL" ~ "B")) %>%
  left_join(letters_groups, by = c("type", "x1_group"))
colors <- rev(brewer.pal(11, "RdBu"))
p <- ggplot(df_long, aes(x = x1_group, y = y, fill = type)) +
  geom_boxplot(
    width = 0.5,
    coef = 0.5,
    position = position_dodge(0.5),
    outlier.shape = NA,
    linewidth = 0.3
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey30", linewidth = 0.3) +
  geom_text(
    data = letters_df,
    aes(x = x1_group, y = 57, label = letter_upper, group = type),
    position = position_dodge(width = 0.5),
    inherit.aes = FALSE, size = 2.5, color = "black"
  ) +
  geom_text(
    data = letters_df,
    aes(x = x1_group, y = 53, label = letter_lower, group = type),
    position = position_dodge(width = 0.5),
    inherit.aes = FALSE, size = 2.3, color = "black"
  ) +
  scale_fill_manual(values = c("Trad" = colors[9], "STL" = colors[3]),
          labels = c("Trad" = "Trad - SSA", "STL" = "STL - SSA")) +
  scale_y_continuous(
  limits = c(-18, 60),
    breaks = seq(0, 40, by = 20),
    expand = c(0, 0)
  ) +
  labs(
    x = "Vegetation type",
    y = expression(Delta ~ tau ~ "(day)"),
    fill = NULL
  ) +
  theme_classic(base_size = 7, base_family = "Arial") +
  theme(
    axis.line = element_blank(),
    axis.ticks = element_line(color = "black", linewidth = 0.4),
    axis.text = element_text(color = "black", size = 7),
    axis.text.x = element_text(margin = margin(t = 2)),
    axis.text.y = element_text(margin = margin(r = 2)),
    axis.title.x = element_text(size = 8, margin = margin(t = 3)),
    axis.title.y = element_text(size = 8, margin = margin(r = 3)),
    legend.position = c(0.25, 0.10),
    legend.key.height = unit(0.3, "cm"),
    legend.key.width = unit(0.2, "cm"),
    legend.background = element_rect(fill = scales::alpha("white", 0.6), color = NA),
    legend.text = element_text(size = 7),
    plot.margin = margin(5, 5, 5, 5),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6)
  )
ggsave("/Volumes/F3c_lulc_gldc.png",
       plot = p, width = 2.2, height = 2.5, dpi = 600,
        bg = "white")