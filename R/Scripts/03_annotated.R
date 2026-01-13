#* 3: Annotated Plots
#+ 3.1: T-tests for all metabolomic features against PGD status
#- 3.1.1: Transform feature table
TFT_annot_transformed <- TFT_annot |> 
  mutate(stage_bin = ifelse(stage_bin == "Early", "Stage I/II", "Stage III/IV"))
#- 3.1.2: Run targeted t-tests
annot_results <- run_targeted_ttests(
  feature_table = TFT_annot_transformed,
  grouping_var = "stage_bin",
  fc_ref_group = "Stage I/II"
)
#- 3.1.3: Add metadata to annotated results
annot_results_w_meta <- annot_results |>
  select(-low_detect_likely) |>
  left_join(TFT_annot_key |> rename(feature = Feature), by = "feature") |>
  arrange(p_value) |>
  select(`Identified Name`, Isomer, log2FC, p_value, p_value_fdr, everything(), -c(unique_vals_no_severe, unique_vals_severe, unique_vals, unique_percentage, n_stage_i_ii, n_stage_iii_iv)) |>
  filter(p_value < 0.05)  
#- 3.1.4: Export as Excel for QC
write.xlsx(annot_results_w_meta, "Outputs/Annotation/annot_results.xlsx")
#! QC done externally and read in as a parameter in the yaml
#+ 3.2: Create Diverging Bar Plots
#- 3.2.1: Prepare data for diverging bars
QC_div <- QC_dedup |>
  filter(diverging_plot == "Y") |>
  select(display_name, log2FC, p_value, main_group) |>
  # Rename groups to match desired labels
  mutate(main_group = case_when(
    main_group == "Steroid Metabolism" ~ "Steroid Flux",
    TRUE ~ main_group
  )) |>
  # Set custom group order from top to bottom
  mutate(main_group = factor(main_group, levels = c(
    "Nucleotide Flux",
    "Epigenetic Signaling",
    "Protein/AA Turnover",
    "Lipid Remodeling",
    "Membrane Integrity",
    "Glycan Defense",
    "Steroid Flux",
    "Immune/Signaling",
    "Fatty Acid Oxidation",
    "Bioenergetic Flux",
    "Redox Homeostasis"
  )))
#- 3.2.2: Create Plot
div_bars <- plot_diverging_bars(QC_div, 
  group_ordering = TRUE, 
  add_group_labels = TRUE,
  max_features = 55,  
  fc_threshold = 0,
  x_max = 5.05,
  lower_expand = 0.0001,
  label_pos = 0.78)
#+ 3.3: Create individual feature plots
#- 3.3.1: Subset to features for scatter plots
QC_scatter <- QC_dedup |>
  filter(scatter_plot == "Y")
#- 3.3.2: Select only scatter plot features
TFT_scatter <- TFT_annot_transformed |>
  select(ID, stage_bin, all_of(QC_scatter$feature))
#- 3.3.3: Create relevant scatter plots
stage_feature_plots <- plot_stage_targeted(
  feature_table = TFT_annot_transformed,  # Use the same transformed data
  metadata_table = QC_scatter,
  include_individual_points = TRUE,
  undo_log = TRUE,
  text_scale = 0.6,
  use_identified_name = TRUE
)
#+ 3.4: Make Correlations Matrix (Full Version)
#- 3.4.1: Subset to selected features for matrix
QC_matrix_full <- QC_dedup |>
  filter(!is.na(feature))
#- 3.4.2: Create data for correlation matrix
TFT_cor_matrix_data_full <- TFT_annot_transformed |>
  select(ID, all_of(QC_matrix_full$feature))
#- 3.4.3: Create correlation matrix plot
corr_mat <- plot_corr_matrix(
  p_threshold = 0.05,
  feature_table = TFT_cor_matrix_data_full,
  metadata_table = QC_matrix_full,
  output_path = "Outputs/Figures/Raw/corr_matrix_full.png",
  show_labels = FALSE
)
#- 3.4.4: Inspect tibble of results
corr_mat$metabolite_stats  |>
  arrange(desc(n_negative_corr), n_positive_corr)
#+ 3.5: Make Correlations Matrix (Curated Version)
#- 3.5.1: Subset to selected features for matrix
QC_matrix <- QC_dedup |>
  filter(matrix_plot == "Y")
#- 3.5.2: Create data for correlation matrix
TFT_cor_matrix_data <- TFT_annot_transformed |>
  select(ID, all_of(QC_matrix$feature))
#- 3.5.3: Create correlation matrix plot
corr_mat <- plot_corr_matrix(
  p_threshold = 0.05,
  feature_table = TFT_cor_matrix_data,
  metadata_table = QC_matrix,
  output_path = "Outputs/Figures/Raw/p3A.png",
  width = 6,
  height = 6
)
