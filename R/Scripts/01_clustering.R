#* 1: PCA and PLS-DA Analysis (Computation Only)
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
#+ 1.3: Heatmaps
{
source("make_heatmap.R")
result <- make_heatmap(
  data = UFT_filtered,
  id_col = "ID",
  feature_selector = "ttest",
  annotation_col = "stage_bin",
  top_features = 200,
  annotation_colors = c("Early" = "#113d6a", "Advanced" = "#800017"
  )
)
print_to_png(result$heatmap_plot, "UFT_stage_heatmap.png", width = 6, height = 4, dpi = 400)
}
