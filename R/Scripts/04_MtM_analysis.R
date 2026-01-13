#* 4: Metabolite to Metabolome Analysis
#+ 4.1: Metabolite to Metabolome Correlation Analysis
#- 4.1.1: Extract feature list
corr_metabolites <- QC_matrix |>
  pull(feature)
#- 4.1.2: Create a lookup: feature -> display_name
feature_lookup <- QC_matrix |>
  select(feature, display_name) |>
  deframe()
#- 4.1.3: Run correlations
pathway_enrichment_correlations <- mummichog_correlations(
  data = UFT_full,
  target_metabolites = corr_metabolites,
  feature_lookup = feature_lookup
)
#+ 4.2: Run Mummichog on Correlation Results
#- 4.2.1: Run Mummichog for each metabolite (MFN)
mummichog_correlations_mfn <- imap(pathway_enrichment_correlations, function(cor_result, metabolite_name) {
  # Create safe directory name from metabolite name
  safe_name <- str_replace_all(metabolite_name, "[^A-Za-z0-9_-]", "_")
  # Run mummichog
  run_mummichog_analysis(
    ttest_results = cor_result,
    output_dir = file.path("Outputs/mummichog/correlations_MFN", safe_name),
    database = "hsa_mfn",
    instrumentOpt = 5.0,
    msModeOpt = "mixed",
    force_primary_ion = "yes",
    p_threshold = 0.05
  )
})
#- 4.2.2: Run Mummichog for each metabolite (KEGG)
mummichog_correlations_kegg <- imap(pathway_enrichment_correlations, function(cor_result, metabolite_name) {
  # Create safe directory name from metabolite name
  safe_name <- str_replace_all(metabolite_name, "[^A-Za-z0-9_-]", "_")
  # Run mummichog
  run_mummichog_analysis(
    ttest_results = cor_result,
    output_dir = file.path("Outputs/mummichog/correlations_KEGG", safe_name),
    database = "hsa_kegg",
    instrumentOpt = 5.0,
    msModeOpt = "mixed",
    force_primary_ion = "yes",
    p_threshold = 0.05
  )
})
#- 4.2.3: Build JSON file paths for all metabolites
correlation_json_paths <- imap(pathway_enrichment_correlations, function(cor_result, metabolite_name) {
  safe_name <- str_replace_all(metabolite_name, "[^A-Za-z0-9_-]", "_")
  list(
    mfn = file.path("Outputs/mummichog/correlations_MFN", safe_name, "scattermum.json"),
    kegg = file.path("Outputs/mummichog/correlations_KEGG", safe_name, "scattermum.json")
  )
})
#+ 4.3: MFN Correlation Enrichment Summary
#- 4.3.1: Import all MFN results
correlation_mfn_tibbles <- imap(correlation_json_paths, function(paths, metabolite_name) {
  read_mummichog_json(paths$mfn) |>
    mutate(metabolite = metabolite_name)
})
#- 4.3.2: Create MFN inspection table
correlation_mfn_inspect <- bind_rows(correlation_mfn_tibbles) |>
  select(metabolite, pathway, enrichment, p_value) |>
  arrange(metabolite, desc(enrichment))
#- 4.3.3: Print MFN summary statistics
cat(
  "\n", strrep("=", 60), "\n",
  "MFN CORRELATION PATHWAY ENRICHMENT SUMMARY\n",
  strrep("=", 60), "\n\n",
  "📊 MFN Database:\n",
  "   Total pathways:   ", nrow(correlation_mfn_inspect), "\n",
  "   Metabolites:      ", length(unique(correlation_mfn_inspect$metabolite)), "\n",
  "   P-value range:    ", min(correlation_mfn_inspect$p_value), " to ", max(correlation_mfn_inspect$p_value), "\n",
  "   Enrichment range: ", round(min(correlation_mfn_inspect$enrichment), 2), " to ", round(max(correlation_mfn_inspect$enrichment), 2), "\n\n",
  strrep("=", 60), "\n\n"
)
#+ 4.4: KEGG Correlation Enrichment Summary
#- 4.4.1: Import all KEGG results
correlation_kegg_tibbles <- imap(correlation_json_paths, function(paths, metabolite_name) {
  read_mummichog_json(paths$kegg) |>
    mutate(metabolite = metabolite_name)
})
#- 4.4.2: Create KEGG inspection table
correlation_kegg_inspect <- bind_rows(correlation_kegg_tibbles) |>
  select(metabolite, pathway, enrichment, p_value) |>
  arrange(metabolite, desc(enrichment))
#- 4.4.3: Print KEGG summary statistics
cat(
  "\n", strrep("=", 60), "\n",
  "KEGG CORRELATION PATHWAY ENRICHMENT SUMMARY\n",
  strrep("=", 60), "\n\n",
  "📊 KEGG Database:\n",
  "   Total pathways:   ", nrow(correlation_kegg_inspect), "\n",
  "   Metabolites:      ", length(unique(correlation_kegg_inspect$metabolite)), "\n",
  "   P-value range:    ", min(correlation_kegg_inspect$p_value), " to ", max(correlation_kegg_inspect$p_value), "\n",
  "   Enrichment range: ", round(min(correlation_kegg_inspect$enrichment), 2), " to ", round(max(correlation_kegg_inspect$enrichment), 2), "\n\n",
  strrep("=", 60), "\n\n"
)
#- 4.4.4: Create combined inspection table
correlation_combined_inspect <- bind_rows(
  correlation_mfn_inspect |> mutate(database = "MFN"),
  correlation_kegg_inspect |> mutate(database = "KEGG")
) |>
  arrange(metabolite, database, desc(enrichment))
#+ 4.5: Visualize MFN Correlation Enrichment
#- 4.5.1: Prepare data for plot_mummichog_columns
correlation_mfn_plot_data <- correlation_mfn_inspect |>
  rename(
    pathway_name = pathway,
    Comparisons = metabolite,
    enrichment_factor = enrichment,
    p_fisher = p_value  # Already -log10(p) from mummichog
  )
#- 4.5.2: Inspect data ranges for plot parameters
cat(
  "\n", strrep("=", 60), "\n",
  "MFN CORRELATION PLOT DATA RANGES\n",
  strrep("=", 60), "\n\n",
  "📊 Plot parameters:\n",
  "   P-value range:    ", min(correlation_mfn_plot_data$p_fisher), " to ", max(correlation_mfn_plot_data$p_fisher), "\n",
  "   Enrichment range: ", round(min(correlation_mfn_plot_data$enrichment_factor), 2), " to ", round(max(correlation_mfn_plot_data$enrichment_factor), 2), "\n",
  "   Pathways:         ", length(unique(correlation_mfn_plot_data$pathway_name)), "\n",
  "   Metabolites:      ", length(unique(correlation_mfn_plot_data$Comparisons)), "\n",
  strrep("=", 60), "\n\n"
)
#- 4.5.3: Create multi-column enrichment plot
correlation_mfn_plot <- plot_mummichog_columns(
  enrichment_data = correlation_mfn_plot_data,
  p_display_threshold = 0.25,
  color_scale = "rb",
  p_color = c(0.25, 0.05, 0.01),
  save_path = "Outputs/Figures/Raw/p3B.png",
  plot_width = 8,
  plot_height = 7,
  detect_thresh = 10,
  detect_sig_thresh = 1,
  hccol = TRUE, 
  hcrow = TRUE,
  dpi = 1000,
  flip_col = TRUE,
  flip_row = FALSE,
  background = "transparent"
)
