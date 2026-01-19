#* 1: Volcano and Heatmaps
#+ 1.1: Volcano Plots on UFT
#- 1.1.1: Run Volcano Analysis
volcano_data <- run_volcano(
  data = UFT_filtered,
  group_var = "stage_bin",
  patient_var = "ID",
  group_levels = c("Early", "Advanced"),
  fc_threshold = log2(1.5),
  p_threshold = 0.05
)
#- 1.1.2: Create Volcano Plot
volcano <- plot_volcano(
  volcano_results = volcano_data,
  x_limits = c(-6, 6),
  y_limits = c(0, 8),
  down_color = "#113d6a"
)
#+ 1.3: Create Heatmaps
#- 1.3.0: Add T_stage_comp to data
UFT_filtered_test <- UFT_filtered |>
  left_join(tumor_pathology |> select(ID, T_stage_comp), by = "ID") |>
  select(ID, T_stage_comp, everything())
#- 1.3.1: Using MAD method
mad_500 <- make_heatmap(
  data = UFT_filtered_test,
  alt_variables = "T_stage_comp",
  stage_colors = c(Early = "#889eb5", Advanced = "#800017"),
  feature_selector = "mad",
  top_features = 500,
  print_preview = FALSE,
  print_scale = FALSE
)
# Visual breakdown of clades...
# Clade 1 = 21 samples, 1/21 advanced; of T category, 6/21 T2, 5/21 T3, remaining T1
# Clade 2 = 34 samples, 5/34 advanced; of T category, 5 T4, 14 T3, 7 T2, remaining T1