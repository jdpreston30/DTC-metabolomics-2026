# Source the pathway name cleaning function
source("R/Utilities/Preprocessing/clean_pathway_names.R")

#+ 3.1: Run t-tests for pairwise variant comparisons
  #- 3.1.1: FV-PTC vs PTC comparison
    variant_assignments_FVPTC_PTC <- UFT_mummichog %>%
      dplyr::filter(Variant != "FTC") %>%
      dplyr::select(Patient_ID, Variant)
    FVPTC_vs_PTC <- mummichog_ttests(
      data = UFT_mummichog %>% dplyr::filter(Variant != "FTC"),
      group_assignments = variant_assignments_FVPTC_PTC,
      group_column = "Variant",
      output_filename = "FVPTC_vs_PTC.csv",
      group1_value = "PTC",
      group2_value = "FV-PTC"
    )
  #- 3.1.2: FTC vs PTC comparison
    variant_assignments_FTC_PTC <- UFT_mummichog %>%
      dplyr::filter(Variant != "FV-PTC") %>%
      dplyr::select(Patient_ID, Variant)
    FTC_vs_PTC <- mummichog_ttests(
      data = UFT_mummichog %>% dplyr::filter(Variant != "FV-PTC"),
      group_assignments = variant_assignments_FTC_PTC,
      group_column = "Variant",
      output_filename = "FTC_vs_PTC.csv",
      group1_value = "PTC",
      group2_value = "FTC"
    )
  #- 3.1.3: FTC vs FV-PTC comparison
    variant_assignments_FTC_FVPTC <- UFT_mummichog %>%
      dplyr::filter(Variant != "PTC") %>%
      dplyr::select(Patient_ID, Variant)
    FTC_vs_FVPTC <- mummichog_ttests(
      data = UFT_mummichog %>% dplyr::filter(Variant != "PTC"),
      group_assignments = variant_assignments_FTC_FVPTC,
      group_column = "Variant",
      output_filename = "FTC_vs_FVPTC.csv",
      group1_value = "FV-PTC",
      group2_value = "FTC"
    )
#+ 3.2: Run t-tests for clusters
  #- 3.2.1: Pull cluster assignments data
    cluster_data <- variance_study %>%
      select(-any_of(feature_cols)) %>%
      left_join(tumor_pathology %>% select(Patient_ID, T_computed), by = "Patient_ID") %>%
      select(Patient_ID, Cluster)
  #- 3.2.2: Add cluster assignments and run ttests
    cluster_assignments <- UFT_mummichog %>%
      left_join(cluster_data, by = "Patient_ID") %>%
      dplyr::select(Patient_ID, Cluster)
    cluster_comparison <- mummichog_ttests(
      data = UFT_mummichog,
      group_assignments = cluster_assignments,
      group_column = "Cluster",
      output_filename = "cluster_v_cluster.csv",
      group1_value = "Cluster 1",
      group2_value = "Cluster 2"
    )
#+ 3.3: Pathway Enrichment Plot (MFN Variants)
  #- 3.6.1: Import results from mummichog
    FVPTC_PTC <- read_xlsx("Outputs/Mummichog Outputs/variant_mummichog.xlsx", sheet = "FVPTC_PTC_MFN") %>%
      mutate(
        Comparisons = "FVPTC_PTC",
        p_fisher = as.numeric(`P(Fisher)`),
        enrichment_factor = Hits.sig / Expected) %>%
        select(Comparisons, pathway_name, p_fisher, enrichment_factor)
    FTC_PTC <- read_xlsx("Outputs/Mummichog Outputs/variant_mummichog.xlsx", sheet = "FTC_PTC_MFN") %>%
      mutate(
        Comparisons = "FTC_PTC",
        p_fisher = as.numeric(`P(Fisher)`),
        enrichment_factor = Hits.sig / Expected) %>%
        select(Comparisons, pathway_name, p_fisher, enrichment_factor)
    FVPTC_FTC <- read_xlsx("Outputs/Mummichog Outputs/variant_mummichog.xlsx", sheet = "FVPTC_FTC_MFN") %>%
      mutate(Comparisons = "FVPTC_FTC",
             p_fisher = as.numeric(`P(Fisher)`),
             enrichment_factor = Hits.sig / Expected) %>%
      select(Comparisons, pathway_name, p_fisher, enrichment_factor)
  #- 3.6.2: Bind rows then filter to important variables
    variant_enrichment <- bind_rows(FVPTC_PTC, FTC_PTC, FVPTC_FTC) %>%
      tidyr::complete(pathway_name, Comparisons) %>%
      filter(p_fisher < 0.05) %>%
      mutate(
        Comparisons = dplyr::case_when(
          Comparisons == "FVPTC_PTC" ~ "PTC vs. FV-PTC",
          Comparisons == "FTC_PTC" ~ "PTC vs. FTC",
          Comparisons == "FVPTC_FTC" ~ "FTC vs. FV-PTC",
          TRUE ~ Comparisons
        ),
        Comparisons = factor(Comparisons, levels = c("PTC vs. FTC", "PTC vs. FV-PTC","FTC vs. FV-PTC")),
        pathway_name = forcats::fct_reorder(pathway_name, enrichment_factor, .fun = max)
      ) %>%
      mutate(enrichment_factor = pmin(enrichment_factor, 5)) %>%
      mutate(pathway_name = clean_pathway_names(pathway_name)) %>%
      mutate(
        pathway_name = factor(
          pathway_name,
          levels = {
            # Get pathways ordered by FTC vs. FV-PTC enrichment factor
            ftc_fvptc_order <- filter(., Comparisons == "FTC vs. FV-PTC" & !is.na(enrichment_factor)) %>%
              arrange(desc(enrichment_factor)) %>%
              pull(pathway_name) %>%
              unique()
            
            # Get any additional pathways not in FTC vs. FV-PTC
            all_pathways <- unique(.$pathway_name)
            remaining_pathways <- setdiff(all_pathways, ftc_fvptc_order)
            
            # Combine: FTC vs. FV-PTC order first, then remaining pathways
            c(ftc_fvptc_order, remaining_pathways)
          }
        )
      )
  #- 3.6.3: Plot
    conflicts_prefer(ggplot2::margin)
    variant_enrichment_plot_MFN <- ggplot(
      variant_enrichment,
      aes(x = 0.5, y = 0.5, size = enrichment_factor, color = p_fisher)
    ) +
      # One dummy row per facet -> avoid the warning
      geom_tile(
        data = data.frame(x = 0.5, y = 0.5),
        aes(x = x, y = y),
        width = 1, height = 1,
        fill = "white", colour = "grey80", linewidth = 0.3,
        inherit.aes = FALSE
      ) +
      geom_point(
        alpha = 0.95, shape = 16, stroke = 0,
        na.rm = TRUE, show.legend = TRUE
      ) +
      facet_grid(
        rows = vars(pathway_name),
        cols = vars(Comparisons),
        switch = "y", drop = FALSE
      ) +
      coord_fixed(clip = "off") +
      scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
      scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +

      # Keep limits ascending; reverse legend order via guide
      scale_size_continuous(
        range = c(5, 10),
        limits = c(0, 5),
        breaks = c(5, 3, 1), # Labels show in this order…
        name = "Enrichment factor",
        guide = guide_legend(reverse = TRUE) # …because we reverse the legend
      ) +

      # Keep p limits ascending; reverse colorbar via guide
      scale_color_gradient(
        low = "#0a2256", high = "#c3dbe9", # Dark (small p) -> light (large p)
        limits = c(0.01, 0.05),
        oob = scales::squish,
        name = "p-value\n",
        guide = guide_colorbar(
          reverse = TRUE, # 0.01 at top, 0.05 at bottom
          barheight = unit(5, "cm"),
          barwidth = unit(0.9, "cm")
        )
      ) +
      labs(x = NULL, y = NULL) +
      theme_minimal(base_family = "Arial") +
      theme(
        text = element_text(family = "Arial"),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.spacing.x = unit(0, "pt"),
        panel.spacing.y = unit(0, "pt"),
        strip.placement = "outside",
        strip.text.x.top = element_text(
          angle = 90, vjust = 0.3, hjust = 0,
          face = "bold", family = "Arial", size = 12, margin = margin(l = -14, b = 5)
        ),
        strip.text.y.left = element_text(
          angle = 0, hjust = 1,
          face = "bold", family = "Arial", size = 12,
          margin = margin(r = 6)
        ),
        legend.title = element_text(size = 12, face = "bold", family = "Arial"),
        legend.text = element_text(size = 12, family = "Arial"),
        plot.margin = margin(t = 20, r = 40, b = 10, l = 40)
      ) +
      coord_cartesian(clip = "off")
  #- 3.6.4: Export Plot
    ggsave(
      "variant_enrichment_plot_MFN.png",
      variant_enrichment_plot_MFN,
      width = length(unique(variant_enrichment$Comparisons)) * 0.3 + 7,
      height = length(unique(variant_enrichment$pathway_name)) * 0.3 + 2,
      units = "in"
    )
#+ 3.3: Pathway Enrichment Plot (KEGG Variants)
  #- 3.6.1: Import results from mummichog
    FVPTC_PTC_KEGG <- read_xlsx("Outputs/Mummichog Outputs/variant_mummichog.xlsx", sheet = "FVPTC_PTC_KEGG") %>%
      mutate(
        Comparisons = "FVPTC_PTC",
        p_fisher = as.numeric(`P(Fisher)`),
        enrichment_factor = Hits.sig / Expected) %>%
        select(Comparisons, pathway_name, p_fisher, enrichment_factor)
    FTC_PTC_KEGG <- read_xlsx("Outputs/Mummichog Outputs/variant_mummichog.xlsx", sheet = "FTC_PTC_KEGG") %>%
      mutate(
        Comparisons = "FTC_PTC",
        p_fisher = as.numeric(`P(Fisher)`),
        enrichment_factor = Hits.sig / Expected) %>%
        select(Comparisons, pathway_name, p_fisher, enrichment_factor)
    FVPTC_FTC_KEGG <- read_xlsx("Outputs/Mummichog Outputs/variant_mummichog.xlsx", sheet = "FVPTC_FTC_KEGG") %>%
      mutate(Comparisons = "FVPTC_FTC",
             p_fisher = as.numeric(`P(Fisher)`),
             enrichment_factor = Hits.sig / Expected) %>%
      select(Comparisons, pathway_name, p_fisher, enrichment_factor)
  #- 3.6.2: Bind rows then filter to important variables
    variant_enrichment_KEGG <- bind_rows(FVPTC_PTC_KEGG, FTC_PTC_KEGG, FVPTC_FTC_KEGG) %>%
      tidyr::complete(pathway_name, Comparisons) %>%
      filter(p_fisher < 0.05) %>%
      mutate(
        Comparisons = dplyr::case_when(
          Comparisons == "FVPTC_PTC" ~ "PTC vs. FV-PTC",
          Comparisons == "FTC_PTC" ~ "PTC vs. FTC",
          Comparisons == "FVPTC_FTC" ~ "FTC vs. FV-PTC",
          TRUE ~ Comparisons
        ),
        Comparisons = factor(Comparisons, levels = c("PTC vs. FTC", "PTC vs. FV-PTC","FTC vs. FV-PTC")),
        pathway_name = forcats::fct_reorder(pathway_name, enrichment_factor, .fun = max)
      ) %>%
      mutate(enrichment_factor = pmin(enrichment_factor, 8)) %>%
      mutate(pathway_name = clean_pathway_names(pathway_name)) %>%
      mutate(
        pathway_name = factor(
          pathway_name,
          levels = {
            # Get pathways ordered by FTC vs. FV-PTC enrichment factor
            ftc_fvptc_order <- filter(., Comparisons == "FTC vs. FV-PTC" & !is.na(enrichment_factor)) %>%
              arrange(desc(enrichment_factor)) %>%
              pull(pathway_name) %>%
              unique()
            
            # Get any additional pathways not in FTC vs. FV-PTC
            all_pathways <- unique(.$pathway_name)
            remaining_pathways <- setdiff(all_pathways, ftc_fvptc_order)
            
            # Combine: FTC vs. FV-PTC order first, then remaining pathways
            c(ftc_fvptc_order, remaining_pathways)
          }
        )
      )
  #- 3.6.3: Plot
    conflicts_prefer(ggplot2::margin)
    variant_enrichment_plot_KEGG <- ggplot(
      variant_enrichment_KEGG,
      aes(x = 0.5, y = 0.5, size = enrichment_factor, color = p_fisher)
    ) +
      # One dummy row per facet -> avoid the warning
      geom_tile(
        data = data.frame(x = 0.5, y = 0.5),
        aes(x = x, y = y),
        width = 1, height = 1,
        fill = "white", colour = "grey80", linewidth = 0.3,
        inherit.aes = FALSE
      ) +
      geom_point(
        alpha = 0.95, shape = 16, stroke = 0,
        na.rm = TRUE, show.legend = TRUE
      ) +
      facet_grid(
        rows = vars(pathway_name),
        cols = vars(Comparisons),
        switch = "y", drop = FALSE
      ) +
      coord_fixed(clip = "off") +
      scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
      scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +

      # Keep limits ascending; reverse legend order via guide
      scale_size_continuous(
        range = c(4, 8),
        limits = c(0, 8),
        breaks = c(8, 6, 4, 2), # Labels show in this order…
        name = "Enrichment factor",
        guide = guide_legend(reverse = TRUE) # …because we reverse the legend
      ) +

      # Keep p limits ascending; reverse colorbar via guide
      scale_color_gradient(
        low = "#0a2256", high = "#c3dbe9", # Dark (small p) -> light (large p)
        limits = c(0.01, 0.05),
        oob = scales::squish,
        name = "p-value\n",
        guide = guide_colorbar(
          reverse = TRUE, # 0.01 at top, 0.05 at bottom
          barheight = unit(5, "cm"),
          barwidth = unit(0.9, "cm")
        )
      ) +
      labs(x = NULL, y = NULL) +
      theme_minimal(base_family = "Arial") +
      theme(
        text = element_text(family = "Arial"),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.spacing.x = unit(0, "pt"),
        panel.spacing.y = unit(0, "pt"),
        strip.placement = "outside",
        strip.text.x.top = element_text(
          angle = 90, vjust = 0.3, hjust = 0,
          face = "bold", family = "Arial", size = 12, margin = margin(l = -14, b = 5)
        ),
        strip.text.y.left = element_text(
          angle = 0, hjust = 1,
          face = "bold", family = "Arial", size = 12,
          margin = margin(r = 6)
        ),
        legend.title = element_text(size = 12, face = "bold", family = "Arial"),
        legend.text = element_text(size = 12, family = "Arial"),
        plot.margin = margin(t = 20, r = 40, b = 10, l = 40)
      ) +
      coord_cartesian(clip = "off")
  #- 3.6.4: Export Plot
    ggsave(
      "variant_enrichment_plot_KEGG.png",
      variant_enrichment_plot_KEGG,
      width = length(unique(variant_enrichment_KEGG$Comparisons)) * 0.3 + 6,
      height = length(unique(variant_enrichment_KEGG$pathway_name)) * 0.3 + 2,
      units = "in"
    )
#+ 3.4: Enrichment Network Plot
  #- Get enrichment factors
      cluster_enrichment_factors <- read_xlsx("Outputs/Mummichog Outputs/variant_mummichog.xlsx", sheet = "FTC_PTC_KEGG") %>%
        mutate(enrichment_factor = Hits.sig / Expected) %>%
        mutate(p_fisher = as.numeric(`P(Fisher)`)) %>%
        select(pathway_name, enrichment_factor)
  #- Import significant FET pathways and 
    KEGG_enrich_raw <- read_xlsx("Outputs/Mummichog Outputs/cluster_mummichog.xlsx", sheet = "EnrichNet_KEGG") %>%
      filter(FET < 0.05) %>%
        left_join(cluster_enrichment_factors, by = "pathway_name") %>%
        mutate(neg_log_p = -log10(FET)) %>%
        select(-FET) %>%
        mutate(
          pathway_name = clean_pathway_names(pathway_name)
        ) %>%
        arrange(enrichment_factor)
  #- Create Plot
    p <- build_enrichment_network(KEGG_enrich_raw, edge_thresh = 0.05, prefer_hsa = TRUE, save_path = "Outputs/Mummichog Outputs/enrichment_network.png", width = 10, height = 8)
