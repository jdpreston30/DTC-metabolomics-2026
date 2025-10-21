#* 3: Annotated Plots
#+ 3.1: T-tests for all metabolomic features against PGD status
#- 3.1.1: Run T-test function for combined TFT
annot_results <- run_targeted_ttests(
  feature_table = TFT_annot,
  grouping_var = "stage_bin",
  fc_ref_group = "Early"
)
#- 3.
annot_results_clean <- annot_results %>%
  filter(low_detect_likely != "Y") %>%
  select(-low_detect_likely) %>%
  left_join(TFT_annot_key %>% select(feature = Feature, `Identified Name`, Isomer), by = "feature") %>%
  arrange(p_value) %>%
  select(`Identified Name`, Isomer, log2FC, p_value, p_value_fdr, everything())
