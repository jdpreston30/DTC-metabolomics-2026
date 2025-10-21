#* 5: Abstract Data Generation
#+ 5.1: Variant Type Distribution
#- 5.1.1: Compute variant counts and percentages
variant_summary <- UFT_filtered %>%
  group_by(Variant) %>%
  summarise(n = n()) %>%
  mutate(percentage = round_half_up(n / sum(n) * 100))

#- 5.1.2: Create variant distribution sentence
variant_counts_sentence <- paste0(
  "Of the ", sum(variant_summary$n), " patients, ",
  paste(
    paste0(variant_summary$n, " (", variant_summary$percentage, "%) had ", variant_summary$Variant),
    collapse = ", "
  ), "."
)

#+ 5.2: Stage Distribution  
#- 5.2.1: Compute stage counts and percentages
stage_summary <- UFT_filtered %>%
  group_by(Stage) %>%
  summarise(n = n()) %>%
  mutate(percentage = round_half_up(n / sum(n) * 100)) %>%
  arrange(Stage)

#- 5.2.2: Create stage distribution sentence
stage_counts_sentence <- paste0(
  "Stage distribution included ",
  paste(
    paste0("Stage ", stage_summary$Stage, " (n=", stage_summary$n, ", ", stage_summary$percentage, "%)"),
    collapse = ", "
  ), "."
)

#+ 5.3: Metabolite Feature Counts
#- 5.3.1: Count total detected features by chromatography method
hilic_count_total <- UFT_full %>% 
  select(starts_with("HILIC")) %>% 
  ncol()

c18_count_total <- UFT_full %>% 
  select(starts_with("C18")) %>% 
  ncol()

total_detected_features <- hilic_count_total + c18_count_total

#- 5.3.2: Count QC-filtered features by chromatography method  
hilic_count_filtered <- UFT_filtered %>% 
  select(starts_with("HILIC")) %>% 
  ncol()

c18_count_filtered <- UFT_filtered %>% 
  select(starts_with("C18")) %>% 
  ncol()

total_filtered_features <- hilic_count_filtered + c18_count_filtered

#- 5.3.3: Create metabolite QC sentence
metabolite_qc_sentence <- paste0(
  "A total of ", total_detected_features, " metabolomic features were initially detected (",
  hilic_count_total, " HILIC and ", c18_count_total, " C18 features). ",
  "After quality control filtering, ", total_filtered_features, " features remained (",
  hilic_count_filtered, " HILIC and ", c18_count_filtered, " C18 features)."
)

#+ 5.4: Volcano Plot Statistics
#- 5.4.1: Extract volcano plot results (assuming you have volcano analysis results)
# Note: Adjust the volcano data object name based on your actual volcano analysis results
if (exists("volcano_results") || exists("pathway_enrichment_ttests")) {
  
  # If using pathway_enrichment_ttests results
  if (exists("pathway_enrichment_ttests")) {
    volcano_data <- pathway_enrichment_ttests$results
    p_threshold <- 0.05
    fc_threshold <- log2(1.5)  # 1.5-fold change
  } else {
    # Adjust these based on your actual volcano results object
    volcano_data <- volcano_results$volcano_data
    p_threshold <- volcano_results$p_threshold
    fc_threshold <- volcano_results$fc_threshold
  }
  
  # Total significantly different features (p < 0.05)
  total_sig <- volcano_data %>%
    filter(p_value < p_threshold) %>%
    nrow()
  
  # Features that are significantly different AND meet fold change threshold (≥1.5-fold)
  sig_up_fc <- volcano_data %>%
    filter(p_value < p_threshold & log2_fc >= fc_threshold) %>%
    nrow()
  
  sig_down_fc <- volcano_data %>%
    filter(p_value < p_threshold & log2_fc <= -fc_threshold) %>%
    nrow()
  
  total_sig_fc <- sig_up_fc + sig_down_fc
  
  #- 5.4.2: Create volcano results sentence
  volcano_results_sentence <- paste0(
    "Differential analysis identified ", total_sig, " significantly different metabolite features (p < ", p_threshold, "), ",
    "of which ", total_sig_fc, " showed ≥1.5-fold change (", sig_up_fc, " increased and ", sig_down_fc, " decreased)."
  )
  
} else {
  volcano_results_sentence <- "Volcano plot analysis results not available - check volcano analysis object names."
}

#+ 5.5: Stage Binning Summary (Early vs Advanced)
#- 5.5.1: Compute early vs advanced stage distribution
if ("stage_bin" %in% colnames(UFT_filtered)) {
  stage_bin_summary <- UFT_filtered %>%
    group_by(stage_bin) %>%
    summarise(n = n()) %>%
    mutate(percentage = round_half_up(n / sum(n) * 100))
  
  early_n <- stage_bin_summary$n[stage_bin_summary$stage_bin == "Early"]
  advanced_n <- stage_bin_summary$n[stage_bin_summary$stage_bin == "Advanced"]
  early_pct <- stage_bin_summary$percentage[stage_bin_summary$stage_bin == "Early"]
  advanced_pct <- stage_bin_summary$percentage[stage_bin_summary$stage_bin == "Advanced"]
  
  #- 5.5.2: Create stage binning sentence
  stage_binning_sentence <- paste0(
    "For comparative analysis, patients were grouped into early-stage (n=", early_n, ", ", early_pct, "%) ",
    "and advanced-stage (n=", advanced_n, ", ", advanced_pct, "%) categories."
  )
} else {
  stage_binning_sentence <- "Stage binning information not available - check stage_bin column."
}

#+ 5.6: Display Manuscript Sentences
cat(
  "\n",
  strrep("=", 70), "\n",
  "MANUSCRIPT SENTENCES FOR THYROID CANCER METABOLOMICS STUDY\n",
  strrep("=", 70), "\n",
  "\n",
  "\033[1;31mVariant Distribution:\033[0m\n",
  "\033[3;31m", variant_counts_sentence, "\033[0m\n",
  "\n",
  "\033[1;32mStage Distribution:\033[0m\n",
  "\033[3;32m", stage_counts_sentence, "\033[0m\n",
  "\n",
  "\033[1;34mStage Grouping (Early vs Advanced):\033[0m\n",
  "\033[3;34m", stage_binning_sentence, "\033[0m\n",
  "\n",
  "\033[1;35mMetabolite Feature Quality Control:\033[0m\n",
  "\033[3;35m", metabolite_qc_sentence, "\033[0m\n",
  "\n",
  "\033[1;36mDifferential Analysis Results:\033[0m\n",
  "\033[3;36m", volcano_results_sentence, "\033[0m\n",
  "\n",
  strrep("=", 70), "\n",
  "\n"
)

#+ 5.7: Summary Statistics Table
#- 5.7.1: Create summary statistics for abstract/manuscript
summary_stats <- data.frame(
  Metric = c(
    "Total Patients",
    "Variant Types",
    "Stage Distribution",
    "Stage Grouping",
    "Total Features Detected",
    "Features After QC",
    "Significantly Different Features",
    "Features with ≥1.5-fold Change"
  ),
  Value = c(
    sum(variant_summary$n),
    paste(variant_summary$Variant, collapse = ", "),
    paste(paste0("Stage ", stage_summary$Stage, " (n=", stage_summary$n, ")"), collapse = "; "),
    if (exists("early_n")) paste0("Early: ", early_n, "; Advanced: ", advanced_n) else "Not available",
    total_detected_features,
    total_filtered_features,
    if (exists("total_sig")) total_sig else "Not calculated",
    if (exists("total_sig_fc")) total_sig_fc else "Not calculated"
  )
)

cat("\033[1;33mSUMMARY STATISTICS TABLE:\033[0m\n")
print(summary_stats, row.names = FALSE)
cat("\n")
