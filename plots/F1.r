cols_method <- c(
  SSA  = "#D55E00",
  STL  = "#0072B2",
  TRAD = "#009E73"
)
method_order <- c("SSA", "STL", "TRAD")
method_labels <- c("SSA" = "SSA", "STL" = "STL", "TRAD" = "Trad")
file_path <- "/Volumes/all_metrics_1000runs.csv"
df_all <- readr::read_csv(file_path, show_col_types = FALSE)
df_long <- df_all %>%
  pivot_longer(
    cols = -id,
    names_to = c("Method", "part1", "part2"),
    names_sep = "_",
    values_to = "value"
  ) %>%
  mutate(
    metric = paste(part1, part2, sep = "_"),
    Method = factor(Method, levels = method_order),
    value = as.numeric(value)
  ) %>%
  dplyr::select(id, Method, metric, value)
create_violin_plot <- function(data,
                               metric_name,
                               y_label,
                               y_limits,
                               y_breaks,
                               stat_labels,
                               trim_q = c(0.01, 0.99)) {
  plot_df <- data %>%
    dplyr::filter(metric == metric_name) %>%
    dplyr::filter(!is.na(value)) %>%
    dplyr::mutate(Method = factor(Method, levels = method_order))
  plot_df_trim <- plot_df %>%
    dplyr::group_by(Method) %>%
    dplyr::mutate(
      q_low  = quantile(value, trim_q[1], na.rm = TRUE),
      q_high = quantile(value, trim_q[2], na.rm = TRUE)
    ) %>%
    dplyr::filter(value >= q_low & value <= q_high) %>%
    dplyr::ungroup() %>%
    dplyr::select(-q_low, -q_high)
  top_y <- y_limits[2] - 0.01 * diff(y_limits)
  text_df <- tibble(
    Method = factor(names(stat_labels), levels = method_order),
    label  = unname(stat_labels),
    y      = top_y
  )
  bottom_y <- y_limits[1] + 0.06 * diff(y_limits)
  median_df <- plot_df_trim %>%
    dplyr::group_by(Method) %>%
    dplyr::summarise(med = median(value, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(
      y = bottom_y,
      label = sprintf("%.3f", med)
    )
  ggplot(plot_df_trim, aes(x = Method, y = value, fill = Method)) +
    geom_violin(
      scale = "width",
      trim = FALSE,
      alpha = 0.9,
      color = "gray20",
      linewidth = 0.35
    ) +
    geom_boxplot(
      width = 0.25,
      outlier.shape = NA,
      fill = "white",
      color = "gray10",
      linewidth = 0.5
    ) +
    stat_summary(
      fun = mean,
      geom = "point",
      shape = 21,
      size = 2.2,
      stroke = 0.35,
      fill = "white",
      color = "black"
    ) +
    geom_text(
      data = text_df,
      aes(x = Method, y = y, label = label, color = Method),
      inherit.aes = FALSE,
      size = 5.5
    ) +
    geom_text(
      data = median_df,
      aes(x = Method, y = y, label = label),
      inherit.aes = FALSE,
      color = "black",
      size = 4.5
    ) +
    coord_cartesian(ylim = y_limits, clip = "off") +
    scale_x_discrete(labels = method_labels) +
    scale_y_continuous(
      breaks = y_breaks,
      expand = expansion(mult = c(0.02, 0.08))
    ) +
    scale_fill_manual(values = cols_method) +
    scale_color_manual(values = cols_method, guide = "none") +
    labs(
      x = NULL,
      y = y_label
    ) +
    theme_prism(base_size = 12, base_family = "Arial") +
    theme(
      legend.position = "none",
      panel.border = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.8),
      axis.ticks = element_line(color = "black", linewidth = 0.7),
      axis.ticks.length = unit(4, "pt"),
      axis.text.x = element_text(size = 14, color = "black", face = "plain", margin = margin(t = 4)),
      axis.text.y = element_text(size = 14, color = "black", face = "plain", margin = margin(r = 4)),
      axis.title.y = element_text(size = 14, color = "black", face = "plain", margin = margin(r = 6)),
      plot.margin = margin(8, 10, 8, 10)
    )
}
p1 <- create_violin_plot(
  data = df_long,
  metric_name = "trend_rmse",
  y_label = "Trend RMSE",
  y_limits = c(0.38, 0.46),
  y_breaks = seq(0.40, 0.46, by = 0.02),
  stat_labels = c("SSA" = "a", "STL" = "b", "TRAD" = "c")
)
p2 <- create_violin_plot(
  data = df_long,
  metric_name = "seas_rmse",
  y_label = "Seasonal cycle RMSE",
  y_limits = c(0.414, 0.422),
  y_breaks = seq(0.416, 0.422, by = 0.002),
  stat_labels = c("SSA" = "a", "STL" = "b", "TRAD" = "c")
)
p3 <- create_violin_plot(
  data = df_long,
  metric_name = "resid_rmse",
  y_label = "Residual RMSE",
  y_limits = c(0, 0.08),
  y_breaks = seq(0.02, 0.08, by = 0.02),
  stat_labels = c("SSA" = "a", "STL" = "b", "TRAD" = "c")
)
p4 <- create_violin_plot(
  data = df_long,
  metric_name = "delta_tau",
  y_label = expression(Delta ~ tau ~ "(day)"),
  y_limits = c(-5, 15),
  y_breaks = seq(0, 15, by = 5),
  stat_labels = c("SSA" = "a", "STL" = "b", "TRAD" = "c")
)
final_plot <- (p1 + p2) / (p3 + p4) +
  plot_annotation(tag_levels = "a") &
  theme(
    plot.tag = element_text(face = "bold", size = 20, family = "Arial")
  )
print(final_plot)
ggsave(
  filename = "/Volumes/plot/F1.png",
  plot = final_plot,
  width = 9,
  height = 8,
  dpi = 300
)
