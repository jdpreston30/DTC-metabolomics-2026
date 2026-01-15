#* 3: Annotated Plots
#+ 3.1: T-tests for all metabolomic features against stage
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
corr_mat_full <- plot_corr_matrix(
  p_threshold = 0.05,
  feature_table = TFT_cor_matrix_data_full,
  metadata_table = QC_matrix_full,
  output_path = "Outputs/Figures/Raw/corr_matrix_full.png",
  matrix_xlsx = "Outputs/Correlation/corr_matrix_full.xlsx",
  show_labels = FALSE
)
#- 3.4.4: Inspect tibble of results
corr_mat_full$metabolite_stats  |>
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
  matrix_xlsx = "Outputs/Correlation/corr_matrix_curated.xlsx",
  width = 6,
  height = 6
)
#+ 3.6: Negative Correlations (p3B)
#- 3.6.1: Acetyl phosphate vs GMP
AP_GMP <- plot_metabolite_correlation(
  y_metabolite = "Acetyl phosphate",
  x_metabolite = "GMP",
  feature_table = TFT_annot_transformed,
  metadata_table = QC_matrix_full,
  minx = 14, maxx = 22, tickx = 2, 
  miny = 15, maxy = 23, ticky = 2,
  text_scale = 0.6,
  stats_annot_pos = "BL"
)
#- 3.6.2: Adrenaline vs 3-Ketosphingosine
Adr_3KS <- plot_metabolite_correlation(
  y_metabolite = "Adrenaline",
  x_metabolite = "3-Ketosphingosine",
  feature_table = TFT_annot_transformed,
  metadata_table = QC_matrix_full,
  minx = 9, maxx = 25, tickx = 4, 
  miny = 14, maxy = 20, ticky = 2,
  text_scale = 0.6,
  stats_annot_pos = "BL"
)
#- 3.6.3: Adrenaline vs Kynurenine
Adr_Kyn <- plot_metabolite_correlation(
  y_metabolite = "Adrenaline",
  x_metabolite = "Kynurenine†",
  feature_table = TFT_annot_transformed,
  metadata_table = QC_matrix_full,
  minx = 11, maxx = 26, tickx = 3, 
  miny = 14, maxy = 20, ticky = 2,
  text_scale = 0.6,
  stats_annot_pos = "BL"
)
#- 3.6.4: Adrenaline vs SAH
Adr_SAH <- plot_metabolite_correlation(
  y_metabolite = "Adrenaline",
  x_metabolite = "SAH",
  feature_table = TFT_annot_transformed,
  metadata_table = QC_matrix_full,
  minx = 14, maxx = 26, tickx = 3, 
  miny = 14, maxy = 20, ticky = 2,
  text_scale = 0.6,
  stats_annot_pos = "BL"
)
#+ 3.7: Positive Correlations (p3D)
#- 3.7.1: GMP vs Ribose 5-phosphate
GMP_R5P <- plot_metabolite_correlation(
  y_metabolite = "GMP",
  x_metabolite = "Ribose 5-phosphate",
  feature_table = TFT_annot_transformed,
  metadata_table = QC_matrix_full,
  minx = 7, maxx = 21, tickx = 2, 
  miny = 13, maxy = 22, ticky = 3,
  text_scale = 0.6
)
#- 3.7.2: 3-Ketosphingosine vs Oleoyl-DHAP
KS3_ODHAP <- plot_metabolite_correlation(
  y_metabolite = "3-Ketosphingosine",
  x_metabolite = "Oleoyl-DHAP",
  feature_table = TFT_annot_transformed,
  metadata_table = QC_matrix_full,
  minx = 11, maxx = 21, tickx = 2, 
  miny = 9, maxy = 25, ticky = 4,
  text_scale = 0.6
)
#- 3.7.3: 1-Methylnicotinamide vs SAH
MNA1_SAH <- plot_metabolite_correlation(
  y_metabolite = "1-Methylnicotinamide",
  x_metabolite = "SAH",
  feature_table = TFT_annot_transformed,
  metadata_table = QC_matrix_full,
  minx = 14, maxx = 26, tickx = 3, 
  miny = 14, maxy = 28, ticky = 2,
  text_scale = 0.6
)
#- 3.7.4: SAH vs PAPS (New Positive)
PAPS_SAH <- plot_metabolite_correlation(
  y_metabolite = "PAPS",
  x_metabolite = "SAH",
  feature_table = TFT_annot_transformed,
  metadata_table = QC_matrix_full,
  minx = 15, maxx = 27, tickx = 3, 
  miny = 4, maxy = 20, ticky = 4,
  text_scale = 0.6
)
#- 3.7.5: Acetylglutamate vs Citrulline
AcGlu_Cit <- plot_metabolite_correlation(
  y_metabolite = "Acetylglutamate",
  x_metabolite = "Citrulline",
  feature_table = TFT_annot_transformed,
  metadata_table = QC_matrix_full,
  minx = 13, maxx = 21, tickx = 2, 
  miny = 15, maxy = 22, ticky = 1,
  text_scale = 0.6
)
#- 3.7.6: Kynurenine vs Serotonin
Kyn_Ser <- plot_metabolite_correlation(
  y_metabolite = "Kynurenine†",
  x_metabolite = "Serotonin",
  feature_table = TFT_annot_transformed,
  metadata_table = QC_matrix_full,
  minx = 11, maxx = 23, tickx = 3, 
  miny = 11, maxy = 26, ticky = 3,
  text_scale = 0.6
)