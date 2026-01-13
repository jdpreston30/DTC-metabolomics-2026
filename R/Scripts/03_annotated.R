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
#+ 3.2: Clean up annotated results
annot_results_clean <- annot_results_w_meta |>
  filter(p_value_fdr < 0.05) |>

# C18_111.020034388403_22.4025506549185
# HILIC_113.034612327566_41.2111075102199
# both of these are uracil, choose the graph that looks better                                           



#+ 3.3: Create individual feature plots
stage_feature_plots <- plot_stage_targeted(
  feature_table = TFT_annot_transformed,  # Use the same transformed data
  metadata_table = annot_results_clean,
  include_individual_points = TRUE,
  undo_log = TRUE,
  text_scale = 0.6,
  use_identified_name = TRUE
)