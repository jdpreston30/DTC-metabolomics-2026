#* 6: Data Not Shown
#+ 6.1: Description of metabolite counts
#- 6.1.1: Count features in full untargeted dataset
{
  hilic_counts_full <- UFT_full |> 
    select(starts_with("HILIC")) |> 
    ncol()
  c18_counts_full <- UFT_full |> 
    select(starts_with("C18")) |> 
    ncol()
  total_untargeted_features <- hilic_counts_full + c18_counts_full
}
#- 6.1.2: Count filtered features 
{
  hilic_count_filtered <- UFT_filtered |> 
    select(starts_with("HILIC")) |> 
    ncol()
  c18_count_filtered <- UFT_filtered |> 
    select(starts_with("C18")) |> 
    ncol()
  total_filtered_features <- hilic_count_filtered + c18_count_filtered
}
#- 6.1.3: Narrative printout
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
#+ 6.2: Volcano Plot Statistics
#- 6.2.1: Extract volcano plot results (assuming you have volcano analysis results)
{
  sig_up_fc <- volcano_data$volcano_data |>
    filter(p_value < 0.05 & log2_fc > log2(1.5)) |>
    nrow()
  sig_down_fc <- volcano_data$volcano_data  |>
    filter(p_value < 0.05 & log2_fc < -log2(1.5)) |>
    nrow()
  
  # FDR-corrected counts
  sig_up_fdr <- volcano_data$volcano_data |>
    filter(p_fdr < 0.05 & log2_fc > log2(1.5)) |>
    nrow()
  sig_down_fdr <- volcano_data$volcano_data  |>
    filter(p_fdr < 0.05 & log2_fc < -log2(1.5)) |>
    nrow()
}
#- 6.2.2: Narrative printout
{
  cat(
    "\n", strrep("=", 60), "\n",
    "VOLCANO PLOT SIGNIFICANT FEATURE COUNTS\n",
    strrep("=", 60), "\n\n",
    "Nominal p-value (p < 0.05, |FC| > 1.5):\n",
    "  Upregulated features:   ", format(sig_up_fc, big.mark = ","), "\n",
    "  Downregulated features: ", format(sig_down_fc, big.mark = ","), "\n",
    "  Total:                  ", format(sig_up_fc + sig_down_fc, big.mark = ","), "\n\n",
    "FDR-corrected (q < 0.05, |FC| > 1.5):\n",
    "  Upregulated features:   ", format(sig_up_fdr, big.mark = ","), "\n",
    "  Downregulated features: ", format(sig_down_fdr, big.mark = ","), "\n",
    "  Total:                  ", format(sig_up_fdr + sig_down_fdr, big.mark = ","), "\n",
    strrep("=", 60), "\n\n"
  )
}
#+ 6.3: Annotated Metabolites Summary
#- 6.3.1: Count annotated metabolites in targeted dataset
{
  hilic_annot <- TFT_annot_transformed |> 
    select(starts_with("HILIC")) |> 
    ncol()
  c18_annot <- TFT_annot_transformed |> 
    select(starts_with("C18")) |> 
    ncol()
  total_annotated <- hilic_annot + c18_annot
}
#- 6.3.2: Count significant features
annot_features_sig <- annot_results_w_meta |>
  select(feature, log2FC, p_value, p_value_fdr) |>
  distinct()
#- 6.3.3: Count significant features by direction
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
#- 6.3.4: Narrative printout
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
#+ 6.4: Figure 2B/C Featured Metabolites
#- 6.4.1: Create tibble of specific features
fig2BC_features <- QC_scatter |>
  select(display_name, log2FC) |>
  arrange(log2FC)
#- 6.4.2: Calculate fold change ranges
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
#- 6.4.3: Narrative printout
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
#+ 6.5: PEA Correlation Summary
# Count significant results per pathway
{
  correlation_mfn_dirs <- list.dirs("Outputs/mummichog/correlations_MFN", recursive = FALSE)
  correlation_mfn_plot_data <- purrr::map_dfr(correlation_mfn_dirs, function(d) {
    csv_path <- file.path(d, "mummichog_pathway_enrichment_mummichog.csv")
    if (!file.exists(csv_path)) return(NULL)
    readr::read_csv(csv_path, show_col_types = FALSE) |>
      dplyr::rename(pathway_name = 1, p_fisher = `P(Fisher)`) |>
      dplyr::mutate(
        p_fisher = -log10(p_fisher),
        metabolite = basename(d)
      )
  })

  correlation_mfn_plot_data |>
    dplyr::group_by(pathway_name) |>
    dplyr::summarise(
      n_significant = sum(p_fisher > -log10(0.05)),
      total_comparisons = dplyr::n()
    ) |>
    dplyr::arrange(desc(n_significant)) |>
    print(n = Inf)
}
#+ 6.6: Abbreviations for figure legends
make_figure_abbrev(abbreviation_tibble)
#+ 6.7: Post-hoc Power Analysis
{
  # Group sample sizes
  group_ns <- UFT_filtered |>
    count(stage_bin)
  n_early    <- group_ns |> filter(stage_bin == "Early")    |> pull(n)
  n_advanced <- group_ns |> filter(stage_bin == "Advanced") |> pull(n)

  # Minimum detectable Cohen's d at 80% power (two-sided alpha = 0.05)
  mde_result <- pwr::pwr.t2n.test(
    n1          = n_early,
    n2          = n_advanced,
    sig.level   = 0.05,
    power       = 0.80,
    alternative = "two.sided"
  )

  # Achieved power at medium effect (Cohen's d = 0.5)
  power_medium <- pwr::pwr.t2n.test(
    n1          = n_early,
    n2          = n_advanced,
    d           = 0.5,
    sig.level   = 0.05,
    alternative = "two.sided"
  )

  # Achieved power at large effect (Cohen's d = 0.8)
  power_large <- pwr::pwr.t2n.test(
    n1          = n_early,
    n2          = n_advanced,
    d           = 0.8,
    sig.level   = 0.05,
    alternative = "two.sided"
  )

  # Observed Cohen's d for FDR-significant annotated features (from raw data)
  fdr_sig_features <- annot_results_w_meta |>
    filter(p_value_fdr < 0.05) |>
    pull(feature)

  observed_cohens_d <- purrr::map_dfr(fdr_sig_features, function(feat) {
    vals <- TFT_annot_transformed |>
      select(stage_bin, value = all_of(feat)) |>
      filter(!is.na(value))
    g1 <- vals |> filter(stage_bin == "Stage I/II")  |> pull(value)
    g2 <- vals |> filter(stage_bin == "Stage III/IV") |> pull(value)
    n1 <- length(g1); n2 <- length(g2)
    if (n1 < 2 || n2 < 2) return(NULL)
    sp <- sqrt(((n1 - 1) * var(g1) + (n2 - 1) * var(g2)) / (n1 + n2 - 2))
    tibble(feature = feat, cohens_d = abs(mean(g2) - mean(g1)) / sp)
  })

  cat(
    "\n", strrep("=", 60), "\n",
    "POST-HOC POWER ANALYSIS (Two-sample Welch t-test)\n",
    strrep("=", 60), "\n\n",
    "Sample sizes:\n",
    "  Early-stage:    n =", n_early, "\n",
    "  Advanced-stage: n =", n_advanced, "\n\n",
    "Minimum detectable effect (80% power, alpha = 0.05):\n",
    "  Cohen's d =", round(mde_result$d, 3), "\n\n",
    "Achieved power:\n",
    "  Medium effect (d = 0.5): ", round(power_medium$power * 100, 1), "%\n",
    "  Large effect  (d = 0.8): ", round(power_large$power  * 100, 1), "%\n\n",
    "Observed Cohen's d for FDR-significant annotated features (n =",
      nrow(observed_cohens_d), "):\n",
    "  Median: ", round(median(observed_cohens_d$cohens_d), 3), "\n",
    "  Min:    ", round(min(observed_cohens_d$cohens_d),    3), "\n",
    "  Max:    ", round(max(observed_cohens_d$cohens_d),    3), "\n",
    strrep("=", 60), "\n\n"
  )

  # Narrative result
  cat(
    "\n", strrep("-", 60), "\n",
    "NARRATIVE RESULT\n",
    strrep("-", 60), "\n\n",
    "No a priori power calculation was performed, as samples were",
    "derived from an existing institutional biobank. Post-hoc power",
    "analysis using the observed group sizes (Stage I/II: n =", n_early,
    "; Stage III/IV: n =", n_advanced, ") indicated",
    round(power_medium$power * 100, 1), "% power to detect a medium",
    "effect (Cohen's d = 0.5) and", round(power_large$power * 100, 1),
    "% power to detect a large effect (Cohen's d = 0.8) at a two-sided",
    "alpha of 0.05. The minimum detectable effect at 80% power was",
    "Cohen's d =", round(mde_result$d, 3), ". Among FDR-significant",
    "annotated features (n =", nrow(observed_cohens_d), "), observed",
    "Cohen's d values ranged from", round(min(observed_cohens_d$cohens_d), 2),
    "to", round(max(observed_cohens_d$cohens_d), 2),
    "(median =", paste0(round(median(observed_cohens_d$cohens_d), 2), "),"),
    "consistent with medium-to-large effect sizes. The limited",
    "advanced-stage sample size reflects the rarity of advanced",
    "differentiated thyroid cancer at our institution; underpowering",
    "in this context is expected to increase the false negative rate",
    "rather than inflate false positives, and all reported findings",
    "are presented as discovery-level observations warranting",
    "prospective validation in larger cohorts.\n\n",
    strrep("-", 60), "\n\n"
  )
}
