#* 5: Abstract Data Generation
#+ 5.1: Variant Type Distribution
#- 5.1.0: Join FT with path data
  path_joined <- UFT_filtered %>%
    left_join(tumor_pathology, by = "ID")
#- 5.1.1: Compute variant counts and percentages
  variant_summary <- path_joined %>%
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
stage_summary <- path_joined %>%
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

#- 5.3.2: Count QC-filtered features by chromatography method  
hilic_count_filtered <- UFT_filtered %>% 
  select(starts_with("HILIC")) %>% 
  ncol()
c18_count_filtered <- UFT_filtered %>% 
  select(starts_with("C18")) %>% 
  ncol()
total_filtered_features <- hilic_count_filtered + c18_count_filtered
#+ 5.4: Volcano Plot Statistics
#- 5.4.1: Extract volcano plot results (assuming you have volcano analysis results)
total_sig <- volcano_data$volcano_data %>%
  filter(p_value < 0.05) %>%
  nrow()
sig_up_fc <- volcano_data$volcano_data %>%
  filter(p_value < p_threshold & log2_fc > 0) %>%
  nrow()
sig_down_fc <- volcano_data$volcano_data  %>%
  filter(p_value < p_threshold & log2_fc < 0) %>%
  nrow()

#+ 5.5: Stage Binning Summary (Early vs Advanced)
#- 5.5.1: Compute early vs advanced stage distribution
mfn_inspect %>%
  filter(p_value >= -log10(0.05))
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
