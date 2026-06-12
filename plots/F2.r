data_path <- "/Volumes/Example1/"
files <- list.files(path = data_path, pattern = "*.csv", full.names = TRUE)
read_single_csv <- function(file_path) {
    df <- read.csv(file_path, header = TRUE)
    file_name <- basename(file_path)
    test_id <- str_extract(file_name, "^\\d{3}")
    type_id <- str_extract(file_name, "(?<=_)\\d{5}")
    ideal_val <- df[4, 1]
    method_vals <- df[1:3, 1] - ideal_val
    df_diff <- data.frame(
        x = method_vals,
        Method = c("ssa", "stl", "trad"),
        Test = as.integer(test_id),
        Type = as.integer(type_id)
    )
    return(df_diff)
}
## example 1,2,3,5,6
all_data <- map_dfr(files, read_single_csv) %>%
  mutate(Method = recode(Method, ssa = "SSA", stl = "STL", trad = "Trad"))
df_median <- all_data %>%
    dplyr::select(Type, Method, x) %>%
    dplyr::group_by(Type, Method) %>%
    dplyr::summarise(Median = median(x, na.rm = TRUE), .groups = "drop")
    fill_colors <- c("SSA" = "#D55E00", "STL" = "#0072B2", "Trad" = "#009E73")
    point_colors <- c("SSA" = "#E69F00", "STL" = "#1F78B4", "Trad" = "#66C2A5")
    point_shapes <- c("SSA" = 21, "STL" = 22, "Trad" = 23)
    ggplot(all_data, aes(x = factor(Type), y = x)) +
      geom_jitter(aes(color = Method), width = 0.2, alpha = 0.5, size = 1.3) +
      geom_point(data = df_median,mapping = aes(x = Type, y = Median, shape = Method, fill = Method),inherit.aes = FALSE,size = 4, stroke = 0.3) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray30", linewidth = 0.6) +
      labs(
        x = "Trend amplitude",
        #x = "Seasonal modulation",
        #x = expression("AR(1) coeff. " * phi),
        #x = "Sampling number per year",
        #x = "Time series length (year)",
        #x = "Random noise SD",
        y = expression(Delta * tau ~ "(day)")
      ) +
      scale_fill_manual(values = fill_colors) +
      scale_color_manual(values = point_colors) +
      scale_shape_manual(values = point_shapes) +
## example 1, 2, 3, 5, 6
       scale_y_continuous(breaks = seq(-5, 15, by = 5),expand = expansion(mult = c(0.1, 0.15)))+
       coord_cartesian(ylim = c(-5, 15))+
## example 4
        #   scale_y_continuous(breaks = seq(-5, 25, by = 5),expand = expansion(mult = c(0.1, 0.15)))+
        #   coord_cartesian(ylim = c(-5, 25))+
      #scale_y_continuous(breaks = c(-5, 0, 5, 10, 90, 100), expand = expansion(mult = c(0.05, 0.05))) +
      #scale_y_break(c(12, 85), scales = 0.6) +
      #coord_cartesian(ylim = c(-5, 105)) +
      # trend example1
      scale_x_discrete(labels = c(0.14, 0.16, 0.18, 0.20, 0.22)) +
      # seasonal example2
      #scale_x_discrete(labels = c(0.05, 0.1, 0.3, 0.6, 0.8)) +
      # anomaly example3
      #scale_x_discrete(labels = c(0.4, 0.5, 0.6, 0.7, 0.8)) +
      ## sampling example4
      ##scale_x_discrete(labels = c(24, 36, 46, 73, 180)) +
      # length example5
      #scale_x_discrete(labels = seq(10, 50, 10)) +
      # noise example6
      #scale_x_discrete(labels = c(0, 0.005, 0.01, 0.02, 0.05)) +
      theme_prism(base_size = 4, base_family = "Arial") +
      theme(
        axis.title.x = element_text(size = 14, face = "plain", margin = margin(t = 3)),
        axis.title.y = element_text(size = 14, face = "plain", margin = margin(r = 1)),
        axis.text.x = element_text(size = 12, color = "black", face = "plain", margin = margin(t = 5)),
        axis.text.y = element_text(size = 12, color = "black", face = "plain", margin = margin(r = 5)),
        legend.position = c(0.22, 0.95),
        legend.justification = c(0, 1),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.text = element_text(size = 13),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(0, 0, 0, 0),
        panel.border = element_blank(),
        axis.line = element_line(color = "black", linewidth = 0.8),
        axis.ticks = element_line(color = "black", linewidth = 0.8),
        axis.ticks.length = unit(4, "pt"),
        panel.grid = element_blank()
      ) +
      guides(
        color = "none",
        fill = guide_legend(title = NULL, direction = "horizontal"),
        shape = guide_legend(title = NULL, direction = "horizontal")
      )
ggsave("/Volumes/F2a_example1_trend.png", width = 4, height = 4, dpi = 300)