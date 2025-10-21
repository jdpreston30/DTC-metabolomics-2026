#* 1: PCA and PLS-DA Analysis (Computation Only)
#+ 1.2: Run PLSDA on UFT data
#+ 1.3: Run Volcano Analysis on UFT data (Computation Only)
volcano_data <- run_volcano(
  data = UFT_filtered,
  group_var = "stage_bin",
  patient_var = "ID",
  group_levels = c("Early", "Advanced"),
  fc_threshold = log2(1.5),
  p_threshold = 0.05
)
#+ 1.4: Create Plots
#- 1.4.6: Volcano
volcano <- plot_volcano(
  volcano_results = volcano_data,
  x_limits = c(-6, 6),
  y_limits = c(0, 8),
  down_color = "#113d6a"
)
