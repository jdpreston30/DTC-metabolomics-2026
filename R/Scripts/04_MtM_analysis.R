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
    ttest_results = list(results = cor_result),
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
    ttest_results = list(results = cor_result),
    output_dir = file.path("Outputs/mummichog/correlations_KEGG", safe_name),
    database = "hsa_kegg",
    instrumentOpt = 5.0,
    msModeOpt = "mixed",
    force_primary_ion = "yes",
    p_threshold = 0.05
  )
})