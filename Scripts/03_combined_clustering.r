#* 3: Combined Clustering analysis
  #+ 3.0: Up front variables
    variant_colors <- c("PTC" = "#DF8D0A", "FV-PTC" = "#23744E", "FTC" = "#194992")
    clade_colors <- c("Clade_1" = "#94001E", "Clade_2" = "#03507D")
  #+ 3.1: Combined PCA analysis based on variants
    pca_variant_results <- make_PCA(
      data = UFT_metaboanalyst_log2,
      output_filename = "Combined_pca_variants.svg",
      ellipse_colors = variant_colors
    )
  #+ 3.2: Heatmap with Combined Untargeted
    #- 3.2.1: Run heatmap analysis
      #_Create heatmap and get results including clade analysis
      heatmap_results <- make_heatmap(
        data = UFT_metaboanalyst_log2,
        variant_colors = variant_colors,
        top_features = 500,
        variant_levels = c("FTC", "FV-PTC", "PTC"),
        n_clades = 2,
        output_filename = "Combined_heatmap.svg"
      )
  #+ 3.3: Combined PCA analysis with hierarchical clusters
    #- 3.3.1: Create modified dataset with clade assignments
      #_Create data with clade assignments instead of variants
      UFT_with_clades <- UFT_metaboanalyst_log2 %>%
        dplyr::left_join(
          heatmap_results$clade_df %>% dplyr::select(Patient_ID, Clade),
          by = "Patient_ID"
        ) %>%
        dplyr::mutate(
          Variant = factor(paste0("Clade_", Clade), levels = c("Clade_1", "Clade_2"))
        ) %>%
        dplyr::select(-Clade) %>%
        dplyr::filter(!is.na(Variant))
    #- 3.3.2: Run PCA with clade colorss
      pca_clades_results <- make_PCA(
        data = UFT_with_clades,
        output_filename = "Combined_pca_clades.svg",
        ellipse_colors = clade_colors
      )
  #+ 3.4: T-tests between hierarchical clusters
    #- 3.4.1: Run t-tests between clades for Mummichog analysis
      ttest_results <- mummichog_ttests(
        data = UFT_metaboanalyst_log2,
        group_assignments = heatmap_results$clade_df,
        group_column = "Clade",
        output_filename = "clade_based_ttest_results.csv",
        group1_value = 1,
        group2_value = 2
      )
  #+ 3.5: Run t-tests for pairwise variant comparisons
    #- 3.5.1: FV-PTC vs PTC comparison
      #_Extract variant assignments (excluding FTC)
      variant_assignments_FVPTC_PTC <- UFT_metaboanalyst_log2 %>%
        dplyr::filter(Variant != "FTC") %>%
        dplyr::select(Patient_ID, Variant)
      FVPTC_vs_PTC <- mummichog_ttests(
        data = UFT_metaboanalyst_log2 %>% dplyr::filter(Variant != "FTC"),
        group_assignments = variant_assignments_FVPTC_PTC,
        group_column = "Variant",
        output_filename = "FVPTC_vs_PTC.csv",
        group1_value = "PTC",
        group2_value = "FV-PTC"
      )
    #- 3.5.2: FTC vs PTC comparison
      #_Extract variant assignments (excluding FV-PTC)
      variant_assignments_FTC_PTC <- UFT_metaboanalyst_log2 %>%
        dplyr::filter(Variant != "FV-PTC") %>%
        dplyr::select(Patient_ID, Variant)
      FTC_vs_PTC <- mummichog_ttests(
        data = UFT_metaboanalyst_log2 %>% dplyr::filter(Variant != "FV-PTC"),
        group_assignments = variant_assignments_FTC_PTC,
        group_column = "Variant",
        output_filename = "FTC_vs_PTC.csv",
        group1_value = "PTC",
        group2_value = "FTC"
      )
    #- 3.5.3: FTC vs FV-PTC comparison
      #_Extract variant assignments (excluding PTC)
      variant_assignments_FTC_FVPTC <- UFT_metaboanalyst_log2 %>%
        dplyr::filter(Variant != "PTC") %>%
        dplyr::select(Patient_ID, Variant)
      FTC_vs_FVPTC <- mummichog_ttests(
        data = UFT_metaboanalyst_log2 %>% dplyr::filter(Variant != "PTC"),
        group_assignments = variant_assignments_FTC_FVPTC,
        group_column = "Variant",
        output_filename = "FTC_vs_FVPTC.csv",
        group1_value = "FV-PTC",
        group2_value = "FTC"
      )
  #+ 3.6: Pathway Enrichment Plot
    #- 3.6.1: Import results from mummichog
      FVPTC_PTC <- read_xlsx("Outputs/Mummichog Outputs/variant_mummichog.xlsx", sheet = "FVPTC_PTC") %>%
        mutate(Comparisons = "FVPTC_PTC")
      FTC_PTC <- read_xlsx("Outputs/Mummichog Outputs/variant_mummichog.xlsx", sheet = "FTC_PTC") %>%
        mutate(Comparisons = "FTC_PTC")
      FVPTC_FTC <- read_xlsx("Outputs/Mummichog Outputs/variant_mummichog.xlsx", sheet = "FVPTC_FTC") %>%
        mutate(Comparisons = "FVPTC_FTC")
    #- 3.6.2: Bind rows then filter to important variables
    pathway_enrich_results <- bind_rows(FVPTC_PTC, FTC_PTC, FVPTC_FTC) %>%
      mutate(enrichment_factor = Hits.sig / Expected) %>%
      select(Comparisons,pathway_name, p_fisher, enrichment_factor) %>%
      arrange(desc(enrichment_factor), p_fisher) %>%
      filter(p_fisher < 0.1)
    #- Create the balloon plot

# Fixed color control on *p_fisher*: min = 0.01 (red), mid = 0.125, max = 0.25 (blue)
pmin <- 0.01
pmid <- 0.05
pmax <- 0.1

df <- pathway_enrich_results %>%
  mutate(
    Comparisons = factor(Comparisons, levels = c("FVPTC_PTC", "FTC_PTC", "FVPTC_FTC")),
    pathway_name = forcats::fct_reorder(pathway_name, enrichment_factor, .fun = max)
  )

ggplot(df, aes(
  x = Comparisons,
  y = pathway_name,
  size = enrichment_factor,
  color = p_fisher
)) +
  geom_point(alpha = 2) +
  scale_size_continuous(range = c(5, 20), name = "Enrichment factor") +
  scale_color_gradient2(
    low = "#f00125", # red at p = 0.01
    mid = "#b54d5d", # midpoint at p = 0.125
    high = "#0678dc", # blue at p = 0.25
    midpoint = pmid,
    limits = c(pmin, pmax),
    oob = scales::squish,
    name = "p (Fisher)"
  ) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_family = "Arial") +
  theme(
    panel.grid.major.y = element_blank(),
    axis.text.y = element_text(size = 11, color = "black"),
    axis.text.x = element_text(size = 11, color = "black"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11)
  )
