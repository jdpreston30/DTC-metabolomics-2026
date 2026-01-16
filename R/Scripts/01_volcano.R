#* 1: Volcano Plots
#+ 1.1: Run Volcano Analysis on UFT data (Computation Only)
volcano_data <- run_volcano(
  data = UFT_filtered,
  group_var = "stage_bin",
  patient_var = "ID",
  group_levels = c("Early", "Advanced"),
  fc_threshold = log2(1.5),
  p_threshold = 0.05
)
#+ 1.2: Create Plots
#- 1.2.6: Volcano
volcano <- plot_volcano(
  volcano_results = volcano_data,
  x_limits = c(-6, 6),
  y_limits = c(0, 8),
  down_color = "#113d6a"
)
#+ 1.3: make_heatmap
u()
ttest_1000 <- make_heatmap(
  data = UFT_filtered,
  feature_selector = "ttest",
  top_features = 1000
)
variance_1000 <- make_heatmap(
  data = UFT_filtered,
  feature_selector = "variance",
  top_features = 1000
)
print_to_png(
  plot = ttest_1000$heatmap_plot,
  filename = "ttest.png",
  dpi = 600,
  width = 6,
  height = 6
)
print_to_png(
  plot = variance_1000$heatmap_plot,
  filename = "variance.png",
  dpi = 600,
  width = 6,
  height = 6
)
