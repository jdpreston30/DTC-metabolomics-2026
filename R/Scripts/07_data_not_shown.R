#* 7: Data Not Shown
#+ 7.1: Description of metabolite counts
#- 7.1.1: Count features in full untargeted dataset
{
  hilic_counts_full <- UFT_full |> 
    select(starts_with("HILIC")) |> 
    ncol()
  c18_counts_full <- UFT_full |> 
    select(starts_with("C18")) |> 
    ncol()
  total_untargeted_features <- hilic_counts_full + c18_counts_full
}
#- 7.1.2: Count filtered features 
{
  hilic_count_filtered <- UFT_filtered |> 
    select(starts_with("HILIC")) |> 
    ncol()
  c18_count_filtered <- UFT_filtered |> 
    select(starts_with("C18")) |> 
    ncol()
  total_filtered_features <- hilic_count_filtered + c18_count_filtered
}
#- 7.1.3: Narrative printout
{
  cat(
    "\n", strrep("=", 60), "\n",
    "FULL UNTARGETED DATASET FEATURE COUNTS\n",
    strrep("=", 60), "\n\n",
    "HILIC (positive mode):  ", format(hilic_counts_full, big.mark = ","), "\n",
    "C18 (negative mode):    ", format(c18_counts_full, big.mark = ","), "\n",
    "Total features:         ", format(total_untargeted_features, big.mark = ","), "\n",
    strrep("=", 60), "\n\n"
  )
  cat(
    "\n", strrep("=", 60), "\n",
    "QC-FILTERED UNTARGETED DATASET FEATURE COUNTS\n",
    strrep("=", 60), "\n\n",
    "HILIC (positive mode):  ", format(hilic_count_filtered, big.mark = ","), "\n",
    "C18 (negative mode):    ", format(c18_count_filtered, big.mark = ","), "\n",
    "Total filtered features:", format(total_filtered_features, big.mark = ","), "\n",
    strrep("=", 60), "\n\n"
  )
}
#+ 7.2: Volcano Plot Statistics
#- 7.2.1: Extract volcano plot results (assuming you have volcano analysis results)
{
  sig_up_fc <- volcano_data$volcano_data |>
    filter(p_value < 0.05 & log2_fc > log2(1.5)) |>
    nrow()
  sig_down_fc <- volcano_data$volcano_data  |>
    filter(p_value < 0.05 & log2_fc < -log2(1.5)) |>
    nrow()
}
#- 7.2.2: Narrative printout
{
  cat(
    "\n", strrep("=", 60), "\n",
    "VOLCANO PLOT SIGNIFICANT FEATURE COUNTS\n",
    strrep("=", 60), "\n\n",
    "Significantly upregulated features (FC > 1.5, p < 0.05):   ", format(sig_up_fc, big.mark = ","), "\n",
    "Significantly downregulated features (FC < -1.5, p < 0.05): ", format(sig_down_fc, big.mark = ","), "\n",
    strrep("=", 60), "\n\n"
  )
}
#+ 7.3: Annotated Metabolites Summary
#- 7.3.1: Count annotated metabolites in targeted dataset
{
  hilic_annot <- TFT_annot_transformed |> 
    select(starts_with("HILIC")) |> 
    ncol()
  c18_annot <- TFT_annot_transformed |> 
    select(starts_with("C18")) |> 
    ncol()
  total_annotated <- hilic_annot + c18_annot
}
#- 7.3.2: Count significant features
annot_features_sig <- annot_results_w_meta |>
  select(feature, log2FC, p_value, p_value_fdr) |>
  distinct()
#- 7.3.3: Count significant features by direction
{
  # Plain p-value < 0.05
  sig_positive_plain <- annot_features_sig |>
    filter(p_value < 0.05 & log2FC > 0) |>
    nrow()
  sig_negative_plain <- annot_features_sig |>
    filter(p_value < 0.05 & log2FC < 0) |>
    nrow()
  sig_total_plain <- sig_positive_plain + sig_negative_plain
  sig_pct_positive <- round(100 * sig_positive_plain / sig_total_plain, 1)
  sig_pct_negative <- round(100 * sig_negative_plain / sig_total_plain, 1)
  # FDR-adjusted p-value < 0.05
  sig_positive_fdr <- annot_features_sig |>
    filter(p_value_fdr < 0.05 & log2FC > 0) |>
    nrow()
  sig_negative_fdr <- annot_features_sig |>
    filter(p_value_fdr < 0.05 & log2FC < 0) |>
    nrow()
  sig_total_fdr <- sig_positive_fdr + sig_negative_fdr
}
#- 7.3.4: Narrative printout
{
  cat(
    "\n", strrep("=", 60), "\n",
    "ANNOTATED METABOLITES IN TARGETED DATASET\n",
    strrep("=", 60), "\n\n",
    "HILIC (positive mode):  ", format(hilic_annot, big.mark = ","), "\n",
    "C18 (negative mode):    ", format(c18_annot, big.mark = ","), "\n",
    "Total annotated metabolites: ", format(total_annotated, big.mark = ","), "\n",
    strrep("=", 60), "\n\n"
  )
  cat(
    "\n", strrep("=", 60), "\n",
    "SIGNIFICANT ANNOTATED FEATURES (EARLY VS ADVANCED)\n",
    strrep("=", 60), "\n\n",
    "Plain p-value < 0.05:\n",
    "  Total significant:            ", format(sig_total_plain, big.mark = ","), "\n",
    "  Increased in advanced stage:  ", format(sig_positive_plain, big.mark = ","), 
    " (", sig_pct_positive, "%)\n",
    "  Decreased in advanced stage:  ", format(sig_negative_plain, big.mark = ","), 
    " (", sig_pct_negative, "%)\n\n",
    "FDR-adjusted p-value < 0.05:\n",
    "  Total significant:            ", format(sig_total_fdr, big.mark = ","), "\n",
    "  Increased in advanced stage:  ", format(sig_positive_fdr, big.mark = ","), "\n",
    "  Decreased in advanced stage:  ", format(sig_negative_fdr, big.mark = ","), "\n\n",
    "NARRATIVE:\n",
    "Among these metabolites, ", sig_total_plain, " were significantly different\n",
    "(", sig_pct_positive, "% positive, ", sig_pct_negative, "% negative) and ", 
    sig_total_fdr, " remained significant\nfollowing FDR correction.\n",
    strrep("=", 60), "\n\n"
  )
}
#+ 7.4: Figure 2B/C Featured Metabolites
#- 7.4.1: Create tibble of specific features
fig2BC_features <- QC_scatter |>
  select(display_name, log2FC) |>
  arrange(log2FC)

#- 7.4.2: Calculate fold change ranges
{
  # Separate increased and decreased metabolites
  increased_metabs <- fig2BC_features |>
    filter(log2FC > 0) |>
    mutate(fold_change = 2^log2FC)
  
  decreased_metabs <- fig2BC_features |>
    filter(log2FC < 0) |>
    mutate(fold_change = 2^abs(log2FC))
  
  # Get ranges
  min_increase <- round_1dp(min(increased_metabs$fold_change))
  max_increase <- round_1dp(max(increased_metabs$fold_change))
  
  min_decrease <- round_1dp(min(decreased_metabs$fold_change))
  max_decrease <- round_1dp(max(decreased_metabs$fold_change))
}
#- 7.4.3: Narrative printout
{
  cat(
    "\n", strrep("=", 60), "\n",
    "FIGURE 2B/C FEATURED METABOLITES\n",
    strrep("=", 60), "\n\n"
  )
  
  cat("Increased metabolites:\n")
  print(increased_metabs, n = Inf)
  cat("\n")
  
  cat("Decreased metabolites:\n")
  print(decreased_metabs, n = Inf)
  cat("\n")
  
  cat(
    "NARRATIVE:\n",
    "Notable among the increased metabolites were GMP, AMP, oleate,\n",
    "γ-linolenate, S-adenosyl-L-homocysteine, and kynurenine,\n",
    "which ranged from ", min_increase, " to ", max_increase, "-fold increased.\n\n",
    "Noteworthy decreased metabolites included 2,3-Dihydroxybenzoate,\n",
    "α-Ketoisocaproate, acetyl phosphate, and adrenaline,\n",
    "which ranged from ", min_decrease, " to ", max_decrease, "-fold decreased.\n",
    strrep("=", 60), "\n\n"
  )
}
#+ 7.5: PEA Correlation Summary
# Count significant results per pathway
{
  correlation_mfn_plot_data |>
    group_by(pathway_name) |>
    summarise(
      n_significant = sum(p_fisher > -log10(0.05)),
      total_comparisons = n()
    ) |>
    arrange(desc(n_significant)) |>
    print(n = Inf)
}
#+ 7.6: Abbreviations for figure legends
make_figure_abbrev(abbreviation_tibble)