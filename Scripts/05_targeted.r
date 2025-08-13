#* 5: Targeted  
#+ 5.1: Merge TFT with path data
  #- 5.1.1: Pull relevant path data
    merged_tumor_data_joiner <- merged_tumor_data %>%
      select(Patient_ID, Clade, T, LD, LVI, ETE)
  #- 5.1.2: Join with TFT_metaboanalyst_log2
    TFT_metaboanalyst_log2_path <- TFT_metaboanalyst_log2 %>%
      left_join(merged_tumor_data_joiner, by = "Patient_ID") %>%
      select(Patient_ID, Clade, T, LD, LVI, Variant,ETE, everything())
#+ 5.2: Statistical comparisons using targeted_metabolite_comparison function
  #- 5.2.0: Load the function
    source("R/Utilities/targeted_comp.R")
  #- 5.2.1: T-test comparison between Clades
    clade_results <- targeted_metabolite_comparison(
      data = TFT_metaboanalyst_log2_path,
      grouping_var = "Clade",
      test_type = "t_test",
      exclude_cols = c("Patient_ID", "T", "LD", "Variant")
    )
  #- 5.1.2: ANOVA comparison between Variants
    variant_results <- targeted_metabolite_comparison(
      data = TFT_metaboanalyst_log2_path,
      grouping_var = "Variant",
      test_type = "anova",
      exclude_cols = c("Patient_ID", "Clade", "T", "LD")
    )
  #- 5.1.3: ANOVA comparison between T stages
    t_stage_results <- targeted_metabolite_comparison(
      data = TFT_metaboanalyst_log2_path,
      grouping_var = "T",
      test_type = "anova",
      exclude_cols = c("Patient_ID", "Clade", "LD", "Variant")
    )

#+ 5.3: Filter metabolites and create top 5 dataset
  #- 5.3.1: Apply isotope filtration to results
    # Filter clade results (remove isotopes)
    clade_results_filtered <- clade_results %>%
      filter(!str_detect(Metabolite, "15N|13C"))
    variant_results_filtered <- variant_results %>%
      filter(!str_detect(Metabolite, "15N|13C"))
  #- 5.3.2: Get top 10 metabolites from clade comparison (by FDR p-value)
    top_10_clade_metabolites <- clade_results_filtered %>%
      filter(!is.na(p_value_fdr)) %>%
      slice_head(n = 10) %>%
      pull(Metabolite)
  #- 5.3.3: Get top 10 metabolites from variant comparison (by FDR p-value)
    top_10_variant_metabolites <- variant_results_filtered %>%
      filter(!is.na(p_value_fdr)) %>%
      slice_head(n = 10) %>%
      pull(Metabolite)
  #- 5.3.4: Create top 10 dataset with metadata
    TFT_top10_dataset_clade <- TFT_metaboanalyst_log2_path %>%
      select(Patient_ID, Clade, T, LD, Variant, all_of(top_10_clade_metabolites))
  #- 5.3.5: Create top 10 dataset with metadata for variant
    TFT_top10_dataset_variant <- TFT_metaboanalyst_log2_path %>%
      select(Patient_ID, Clade, T, LD, Variant, all_of(top_10_variant_metabolites))
  #- 5.3.6: Export top 10 datasets
    write.csv(TFT_top10_dataset_clade, "Outputs/TFT_top10_dataset_clade.csv", row.names = FALSE)
    write.csv(TFT_top10_dataset_variant, "Outputs/TFT_top10_dataset_variant.csv", row.names = FALSE)

# Sample sizes for T × Variant × Clade subclasses
n_per_clade_T_variant <- TFT_metaboanalyst_log2_path %>%
  filter(!is.na(T), !is.na(Variant), !is.na(Clade)) %>%
  count(T, Variant, Clade, name = "n") %>%
  arrange(Clade, Variant, T)

write.csv(n_per_clade_T_variant, "Outputs/n_per_clade_T_variant.csv", row.names = FALSE)

n_per_LVI_clade_variant <- TFT_metaboanalyst_log2_path %>%
  filter(!is.na(T), !is.na(Variant), !is.na(Clade)) %>%
  count(LVI, Variant, Clade, name = "n") %>%
  arrange(Clade, Variant, LVI)

ETE_clade_var <- TFT_metaboanalyst_log2_path %>%
  filter(!is.na(T), !is.na(Variant), !is.na(Clade)) %>%
  count(ETE, Variant, Clade, name = "n") %>%
  arrange(Clade, Variant, ETE)
