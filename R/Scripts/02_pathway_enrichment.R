#* 2: Pathway Enrichment Analysis
#+ 2.1: Run Mummichog Analysis
#- 2.1.1: Run Mummichog Ttest function
pathway_enrichment_ttests <- mummichog_ttests(
  data = UFT_full,
  group_column = "stage_bin",
  group1_value = "Early",
  group2_value = "Advanced"
)
#- 2.1.2: Run Mummichog (MFN Only)
mummichog_mfn <- run_mummichog_analysis(
  ttest_results = pathway_enrichment_ttests$results,
  output_dir = "Outputs/mummichog/MFN",
  database = "hsa_mfn",
  instrumentOpt = 5.0,
  msModeOpt = "mixed",
  force_primary_ion = "yes"
)
#- 2.1.3: Run Mummichog (KEGG)
mummichog_kegg <- run_mummichog_analysis(
  ttest_results = pathway_enrichment_ttests$results,
  output_dir = "Outputs/mummichog/KEGG",
  database = "hsa_kegg",
  instrumentOpt = 5.0,
  msModeOpt = "mixed",
  force_primary_ion = "yes"
)
#+ 2.3: Create Pathway Enrichment Plots
#- 2.3.1: Define JSON file paths once
mfn_json_files <- list(
  eva = "Outputs/mummichog/MFN/scattermum.json"
)
kegg_json_files <- list(
  eva = "Outputs/mummichog/KEGG/scattermum.json"
)
combined_json_files <- list(
  mfn = mfn_json_files,
  kegg = kegg_json_files
)
#- 2.3.2: Import tibbles for inspection
mfn_tibbles <- map(mfn_json_files, read_mummichog_json)
kegg_tibbles <- map(kegg_json_files, read_mummichog_json)
combined_tibbles <- list(
  MFN = mfn_tibbles,
  KEGG = kegg_tibbles
)
#- 2.3.3: Make MFN inspection tibble for visualization
mfn_inspect <- mfn_tibbles$eva %>%
  arrange(desc(enrichment)) %>%
  filter(p_value >= 1)
#- 2.3.4: Inspect combined for visualization
combined_inspection <- bind_rows(
  mfn_tibbles$eva %>% mutate(database = "MFN"),
  kegg_tibbles$eva %>% mutate(database = "KEGG")
) %>%
  arrange(enrichment) %>%
  filter(p_value >= 1)
#- 2.3.5: Print min and max of MFN and combined p values and enrichment for visualization scaling
cat(
  "\n", strrep("=", 60), "\n",
  "PATHWAY ENRICHMENT RANGES FOR VISUALIZATION SCALING\n",
  strrep("=", 60), "\n\n",
  "📊 MFN Database:\n",
  "   P-value range:    ", min(mfn_inspect$p_value), " to ", max(mfn_inspect$p_value), "\n",
  "   Enrichment range: ", min(mfn_inspect$enrichment), " to ", max(mfn_inspect$enrichment), "\n\n",
  "📊 Combined (MFN + KEGG):\n",
  "   P-value range:    ", min(combined_inspection$p_value), " to ", max(combined_inspection$p_value), "\n",
  "   Enrichment range: ", min(combined_inspection$enrichment), " to ", max(combined_inspection$enrichment), "\n\n",
  strrep("=", 60), "\n\n"
)
#- 2.3.6: Make MFN only plot
mfn_enrichment_plot <- plot_mummichog_enrichment(
  json_files = mfn_json_files,
  combine_databases = FALSE,
  p_threshold = 0.1,
  enrichment_cap = 2.5,
  size_range = c(.11, 12),
  size_breaks = c(2.5, 2.0, 1.5, 1.0),
  show_legend = TRUE,
  save_path = "Outputs/Figures/Raw/mfn_enrich.png",
  plot_width = 5.9,
  plot_height = 5.4, # ! 6.3 for 3b version
  dpi = 600,
  color_scale = "rb"
)
#+ 2.4: Run Biological Network Analysis
mfn_network <- create_biological_network(
  pathway_csv = "Outputs/mummichog/MFN/mummichog_pathway_enrichment_mummichog.csv",
  min_shared_compounds = 1,
  p_threshold = 0.1,
  max_pathways = 20,
  network_name = "mfn_biological"
)
#+ 2.5: Plot Biological Networks
network_plot <- plot_biological_network(
  network_data = mfn_network,
  output_file = "Outputs/Figures/Raw/mfn_network.png",
  node_size_range = c(4, 14),
  text_size = 4.5,
  show_legend = FALSE,
  plot_width = 10.5,
  plot_height = 10.5,
  dpi = 600,
  seed = 2025,
  variable_edge_thickness = TRUE,
  edge_thickness_range = c(0.3, 3),
  max_distance_from_center = 1.5,
  label_position = "above",
  show_node_numbers = FALSE,
  labels_below = c(3,7,9),
  nudge_labels_vert = list(p1 = 1, p5 = 0, p10 = 0),
  nudge_labels_horiz = list(p1 = 1.1, p5 = -0.2, p7 = -0.05, p9 = 0.5, p10 = 1.1),
  color_scale = "rb"
)
