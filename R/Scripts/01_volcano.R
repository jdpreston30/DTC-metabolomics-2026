#* 1: Volcano and Heatmaps
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
#- 1.3.1: Using MAD method
mad_500 <- make_heatmap(
  data = UFT_filtered,
  feature_selector = "mad",
  top_features = 500,
  print_preview = TRUE,
  print_scale = TRUE
)
#- 1.3.1: Using ttest method
ttest_500 <- make_heatmap(
  data = UFT_filtered,
  feature_selector = "ttest",
  top_features = 500,
  print_preview = TRUE,
  print_scale = TRUE
)