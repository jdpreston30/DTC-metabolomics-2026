#* 9: Comprehensive Results Printout
#+ 9.1: Volcano Plot Results
cat("\n", strrep("=", 80), "\n")
cat("SECTION 1: GLOBAL METABOLIC PROFILING - VOLCANO PLOT ANALYSIS\n")
cat(strrep("=", 80), "\n\n")

{
  # Extract volcano statistics
  volcano_sig <- volcano_data$volcano_data |>
    filter(significant == TRUE)
  
  volcano_up <- volcano_sig |> filter(log2_fc > 0)
  volcano_down <- volcano_sig |> filter(log2_fc < 0)
  
  total_features <- nrow(volcano_data$volcano_data)
  sig_features <- nrow(volcano_sig)
  sig_pct <- round(100 * sig_features / total_features, 1)
  
  cat("Total metabolomic features analyzed:", total_features, "\n")
  cat("Significantly different features (p<0.05, |FC|>1.5):", sig_features, 
      sprintf("(%s%%)\n", sig_pct))
  cat("  - Upregulated in advanced stage:", nrow(volcano_up), "\n")
  cat("  - Downregulated in advanced stage:", nrow(volcano_down), "\n\n")
  
  cat("Top 10 Most Upregulated Features:\n")
  cat(strrep("-", 80), "\n")
  top_up <- volcano_up |> 
    arrange(desc(log2_fc)) |>
    head(10) |>
    mutate(fold_change = round(2^log2_fc, 2))
  print(top_up |> select(feature, log2_fc, fold_change, p_value), n = Inf)
  
  cat("\n\nTop 10 Most Downregulated Features:\n")
  cat(strrep("-", 80), "\n")
  top_down <- volcano_down |>
    arrange(log2_fc) |>
    head(10) |>
    mutate(fold_change = round(2^abs(log2_fc), 2))
  print(top_down |> select(feature, log2_fc, fold_change, p_value), n = Inf)
}

#+ 9.2: Pathway Enrichment Results (MFN)
cat("\n\n", strrep("=", 80), "\n")
cat("SECTION 2: PATHWAY ENRICHMENT ANALYSIS - MFN DATABASE\n")
cat(strrep("=", 80), "\n\n")

{
  mfn_results <- mfn_tibbles$eva |>
    filter(p_value < 0.1) |>
    arrange(p_value)
  
  cat("Number of significantly enriched pathways (p<0.1):", nrow(mfn_results), "\n\n")
  
  cat("Significantly Enriched Pathways (MFN):\n")
  cat(strrep("-", 80), "\n")
  print(mfn_results |> 
    select(pathway_name, enrichment, p_value, overlap_size, pathway_size), 
    n = Inf)
}

#+ 9.3: Pathway Enrichment Results (KEGG)
cat("\n\n", strrep("=", 80), "\n")
cat("SECTION 3: PATHWAY ENRICHMENT ANALYSIS - KEGG DATABASE\n")
cat(strrep("=", 80), "\n\n")

{
  kegg_results <- kegg_tibbles$eva |>
    filter(p_value < 0.1) |>
    arrange(p_value)
  
  cat("Number of significantly enriched pathways (p<0.1):", nrow(kegg_results), "\n\n")
  
  cat("Significantly Enriched Pathways (KEGG):\n")
  cat(strrep("-", 80), "\n")
  print(kegg_results |> 
    select(pathway_name, enrichment, p_value, overlap_size, pathway_size), 
    n = Inf)
}

#+ 9.4: Annotated Metabolites Analysis
cat("\n\n", strrep("=", 80), "\n")
cat("SECTION 4: ANNOTATED METABOLITES - DIFFERENTIAL ABUNDANCE\n")
cat(strrep("=", 80), "\n\n")

{
  # Get statistics
  annot_sig_plain <- annot_features_sig |> filter(p_value < 0.05)
  annot_sig_fdr <- annot_features_sig |> filter(p_value_fdr < 0.05)
  
  cat("Total annotated metabolites tested:", nrow(annot_features_sig), "\n")
  cat("Significant (p<0.05):", nrow(annot_sig_plain), "\n")
  cat("Significant (FDR<0.05):", nrow(annot_sig_fdr), "\n\n")
  
  # Add display names
  annot_with_names <- annot_results_w_meta |>
    filter(p_value < 0.05) |>
    mutate(fold_change = round(2^abs(log2FC), 2),
           direction = ifelse(log2FC > 0, "UP", "DOWN")) |>
    arrange(p_value)
  
  cat("All Significantly Different Annotated Metabolites (p<0.05):\n")
  cat(strrep("-", 80), "\n")
  print(annot_with_names |> 
    select(`Identified Name`, direction, log2FC, fold_change, p_value, p_value_fdr), 
    n = Inf)
}

#+ 9.5: Correlation Matrix Summary
cat("\n\n", strrep("=", 80), "\n")
cat("SECTION 5: METABOLITE CORRELATION MATRIX ANALYSIS\n")
cat(strrep("=", 80), "\n\n")

{
  # Get correlation matrix statistics
  corr_stats <- corr_mat_full$metabolite_stats |>
    arrange(desc(n_positive_corr))
  
  cat("Correlation matrix dimensions:", nrow(QC_matrix_full), "x", nrow(QC_matrix_full), "\n")
  cat("Significance threshold: p<0.05\n\n")
  
  cat("Top Metabolites by Number of Positive Correlations:\n")
  cat(strrep("-", 80), "\n")
  print(corr_stats |> head(20), n = Inf)
  
  cat("\n\nTop Metabolites by Number of Negative Correlations:\n")
  cat(strrep("-", 80), "\n")
  print(corr_stats |> arrange(desc(n_negative_corr)) |> head(20), n = Inf)
}

#+ 9.6: Specific Correlation Pairs
cat("\n\n", strrep("=", 80), "\n")
cat("SECTION 6: KEY METABOLITE CORRELATION PAIRS\n")
cat(strrep("=", 80), "\n\n")

{
  # Define correlation pairs from Figure 3
  cat("NEGATIVE CORRELATIONS (Efficient Metabolism vs. Tumor Growth):\n")
  cat(strrep("-", 80), "\n\n")
  
  neg_pairs <- list(
    c("Acetyl phosphate", "GMP"),
    c("Adrenaline", "3-Ketosphingosine"),
    c("Adrenaline", "Kynurenine"),
    c("Adrenaline", "SAH")
  )
  
  for (pair in neg_pairs) {
    cat(sprintf("%s vs %s:\n", pair[1], pair[2]))
    # Would need actual correlation values - placeholder
    cat("  [Correlation statistics would be calculated here]\n\n")
  }
  
  cat("\nPOSITIVE CORRELATIONS (Growth & Proliferation Metabolites):\n")
  cat(strrep("-", 80), "\n\n")
  
  pos_pairs <- list(
    c("GMP", "Ribose 5-phosphate"),
    c("3-Ketosphingosine", "Oleoyl-DHAP"),
    c("1-Methylnicotinamide", "SAH"),
    c("PAPS", "SAH"),
    c("Acetylglutamate", "Citrulline"),
    c("Kynurenine", "Serotonin")
  )
  
  for (pair in pos_pairs) {
    cat(sprintf("%s vs %s:\n", pair[1], pair[2]))
    cat("  [Correlation statistics would be calculated here]\n\n")
  }
}

#+ 9.7: Metabolite-to-Metabolome Enrichment
cat("\n\n", strrep("=", 80), "\n")
cat("SECTION 7: METABOLITE-TO-METABOLOME ENRICHMENT ANALYSIS\n")
cat(strrep("=", 80), "\n\n")

{
  cat("Key metabolites analyzed for pathway associations:\n\n")
  
  # List the metabolites that were analyzed
  key_metabs <- c("α-Ketoisocaproate", "γ-Linolenate", "11-cis-Retinol",
                  "2,3-Dihydroxybenzoate", "7-DHC", "Acetyl phosphate",
                  "Acetylglutamate", "Adrenaline", "AMP", "Butanoylcarnitine",
                  "GMP", "Histidine", "Kynurenine", "LTD4", "Lysine",
                  "Oleate", "PAPS", "Ribose 5-phosphate", "SAH", "Serotonin")
  
  for (metab in key_metabs) {
    cat(sprintf("- %s\n", metab))
  }
  
  cat("\n[Detailed pathway associations for each metabolite would be listed here]\n")
  cat("[This includes correlations with pathways like:\n")
  cat("  - Fatty acid metabolism\n")
  cat("  - Amino acid metabolism (Ala, Asp, Gly, Ser, Thr)\n")
  cat("  - Nucleotide metabolism\n")
  cat("  - Lipid remodeling\n")
  cat("  - etc.]\n")
}

#+ 9.8: Featured Metabolites Summary
cat("\n\n", strrep("=", 80), "\n")
cat("SECTION 8: FEATURED METABOLITES IN FIGURES\n")
cat(strrep("=", 80), "\n\n")

{
  cat("FIGURE 2B - Increased Metabolites:\n")
  cat(strrep("-", 80), "\n")
  print(increased_metabs, n = Inf)
  
  cat("\n\nFIGURE 2C - Decreased Metabolites:\n")
  cat(strrep("-", 80), "\n")
  print(decreased_metabs, n = Inf)
  
  cat("\n\nFold Change Summary:\n")
  cat(sprintf("Increased metabolites ranged from %s to %s-fold higher in advanced stage\n",
              min_increase, max_increase))
  cat(sprintf("Decreased metabolites ranged from %s to %s-fold lower in advanced stage\n",
              min_decrease, max_decrease))
}

#+ 9.9: Summary Statistics
cat("\n\n", strrep("=", 80), "\n")
cat("SECTION 9: OVERALL SUMMARY STATISTICS\n")
cat(strrep("=", 80), "\n\n")

{
  cat("DATASET COMPOSITION:\n")
  cat("  Total samples:", nrow(UFT_filtered), "\n")
  cat("  Early stage:", sum(UFT_filtered$stage_bin == "Early"), "\n")
  cat("  Advanced stage:", sum(UFT_filtered$stage_bin == "Advanced"), "\n\n")
  
  cat("FEATURE COUNTS:\n")
  cat("  Total untargeted features (pre-QC):", total_untargeted_features, "\n")
  cat("  QC-filtered features:", total_filtered_features, "\n")
  cat("  Annotated metabolites:", total_annotated, "\n\n")
  
  cat("SIGNIFICANCE SUMMARY:\n")
  cat("  Volcano plot significant features:", sig_features, sprintf("(%s%%)\n", sig_pct))
  cat("  Annotated metabolites (p<0.05):", nrow(annot_sig_plain), "\n")
  cat("  Annotated metabolites (FDR<0.05):", nrow(annot_sig_fdr), "\n")
  cat("  MFN enriched pathways (p<0.1):", nrow(mfn_results), "\n")
  cat("  KEGG enriched pathways (p<0.1):", nrow(kegg_results), "\n")
}

cat("\n\n", strrep("=", 80), "\n")
cat("END OF RESULTS PRINTOUT\n")
cat(strrep("=", 80), "\n\n")
