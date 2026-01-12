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
#- 3.1.3: Clean up annotated results
annot_results_clean <- annot_results |>
  select(-low_detect_likely) |>
  left_join(TFT_annot_key |> select(feature = Feature, `Identified Name`, Isomer), by = "feature") |>
  arrange(p_value) |>
  select(`Identified Name`, Isomer, log2FC, p_value, p_value_fdr, everything()) |>
  filter(p_value_fdr < 0.05) |>
  filter(`Identified Name` %in% c(
    "GMP",
    "AMP",
    "(6Z,9Z,12Z)-Octadecatrienoic acid",
    "S-Adenosyl-L-homocysteine"
  )) |>
  mutate(`Identified Name` = ifelse(`Identified Name` == "(6Z,9Z,12Z)-Octadecatrienoic acid", "γ-Linolenic Acid", `Identified Name`)) |>
  mutate(`Identified Name` = ifelse(`Identified Name` == "S-Adenosyl-L-homocysteine", "SAH", `Identified Name`))
#+ 3.2: Source the plotting function
source("R/Utilities/Visualization/plot_stage_targeted.R")
#+ 3.3: Create individual feature plots
stage_feature_plots <- plot_stage_targeted(
  feature_table = TFT_annot_transformed,  # Use the same transformed data
  metadata_table = annot_results_clean,
  include_individual_points = TRUE,
  undo_log = TRUE,
  text_scale = 0.6,
  use_identified_name = TRUE
)