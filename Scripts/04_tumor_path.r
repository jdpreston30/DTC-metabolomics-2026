#* 4: Tumor Pathology Analysis
#! Need to look close at P12 and P13 to see if it is correct to grade them as 4 based on ETE as others with minimal ETE are graded as 3
#+ 4.1: Tumor Path based differences by clusters
  #- 4.1.1: Examine clade based differences
  path_comparison_clades <- ternG(
    data = merged_tumor_data %>% select(-c(Patient_ID,Patient)),
    group_var = "Clade",
    force_ordinal = "T",
    descriptive = FALSE,
    OR_col = FALSE,
    OR_method = "dynamic",
    consider_normality = FALSE,
    print_normality = FALSE
  )
  #- 4.1.2: Examine potential confounding between clade and variant
  summary(lm(LD ~ Clade + Variant_Name, data = merged_tumor_data))
  summary(polr(as.factor(T) ~ Clade + Variant_Name, data = merged_tumor_data, Hess = TRUE))
  #- 4.1.3: Compute T percentages for parts of whole graph
    stage_pct <- merged_tumor_data_export %>%
      count(Clade, T) %>%
      group_by(Clade) %>%
      mutate(
        pct = 100 * n / sum(n)
      ) %>%
      ungroup()
  #- 4.1.3: Export results for prism
    output_csv(stage_pct, "tumor_pathology_clades_stage_pct.csv")
    output_csv(merged_tumor_data %>% select(Clade, LD, T), "tumor_pathology_clades_results.csv")

