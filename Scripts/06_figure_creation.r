#* 6: Figure Creation
  #+ 6.1: Figure 1 - Clustering Analysis Layout
    #- 6.1.1: Fix the border widths
      permanova_1B <- permanova_1B + theme_pub_simple(border_linewidth = 0.5)
      pca_1C <- pca_1C + theme_pub_pca(border_linewidth = 0.5)
      pca_1D <- pca_1D + theme_pub_pca(border_linewidth = 0.5)
      pca_1E <- pca_1E + theme_pub_pca(border_linewidth = 0.5)
      pca_1F <- pca_1F + theme_pub_pca(border_linewidth = 0.5)
    #- 6.1.4: Assemble Figure 1 (Heatmap top 50%, vertical PERMANOVA + 4 PCAs bottom 50%)
      Figure_1 <- patchwork::wrap_plots(
        # Top 50%: Heatmap spanning full width
        heatmap_1A,
        # Bottom 50%: Vertical PERMANOVA (left) + 4 PCAs (2x2 grid, right)
        permanova_1B, pca_1C, pca_1D, pca_1E, pca_1F,
        design = "AAAA\nAAAA\nBBCD\nBBEF",  # Heatmap top 50%, vertical bar + 2x2 PCAs bottom 50%
        heights = c(1, 1, 0.5, 0.5)  # 50% heatmap, 50% bottom split
      ) +
        patchwork::plot_annotation(
          title = "Figure 1",
          tag_levels = "A",
          theme = ggplot2::theme(
            plot.title.position = "plot",
            plot.title = ggplot2::element_text(hjust = 0, face = "bold", family = "Arial", size = 16),
            plot.margin = grid::unit(c(0.3, 0.5, 0.3, 0.5), "in")
          )
        ) &
        ggplot2::theme(
          plot.tag.position = c(0, 0.98),
          plot.tag = ggplot2::element_text(size = 14, face = "bold", vjust = 0, hjust = 0, family = "Arial", color = "black")
        )
    #- 6.1.6: Save Figure 1 to PDF and PNG
      print_to_png(Figure_1, "Figure_1_preview")
  #+ 6.3: Figure 3 - Targeted Metabolomics
    #- 6.3.1 Assemble Figure 3
    Figure_3 <- patchwork::wrap_plots(
      facet_3A, clusters_targ_3B,                # Top 25%: A, B
      variant_targ_3C,
      design = "A\nB\nC",
      heights = c(0.25, 0.25, 0.25, 0.25)
    ) +
      patchwork::plot_annotation(
        title = "Figure 3",
        tag_levels = "A",
        theme = ggplot2::theme(
          plot.title.position = "plot",
          plot.title = ggplot2::element_text(hjust = 0, face = "bold", family = "Arial", size = 16),
          plot.margin = grid::unit(c(0.1, 1.9, 0.1, 1.9), "in")
        )
      ) &
      ggplot2::theme(
        plot.tag.position = c(0, 0.98),
        plot.tag = ggplot2::element_text(size = 14, face = "bold", vjust = 0, hjust = 0, family = "Arial", color = "black")
      )
    #- 6.3.2 Save Figure 3 to PNG
      print_to_png(Figure_3, "Figure_3_preview")
  #+ 6.4: Supplemental Figure 1 - MFN Pathway Enrichment
    #- 6.4.1: Calculate MFN plot dimensions and scale to 8.5 x 11
      # Original proportions that work well
      mfn_width_orig <- length(unique(variant_enrichment$Comparisons)) * 0.3 + 7
      mfn_height_orig <- length(unique(variant_enrichment$pathway_name)) * 0.3 + 2
      
      # Calculate aspect ratio
      aspect_ratio <- mfn_width_orig / mfn_height_orig
      
      # Scale to fit within 8.5 x 11 (leaving margins for title)
      max_width <- 7.5   # 8.5 - 1 for margins
      max_height <- 10   # 11 - 1 for title and margins
      
      # Scale while preserving aspect ratio
      if (aspect_ratio > max_width/max_height) {
        # Width-constrained
        scaled_width <- max_width
        scaled_height <- max_width / aspect_ratio
      } else {
        # Height-constrained  
        scaled_height <- max_height
        scaled_width <- max_height * aspect_ratio
      }
      
    #- 6.4.2: Create Supplemental Figure 1
      Supplemental_Figure_1 <- variant_enrichment_plot_MFN +
        patchwork::plot_annotation(
          title = "Supplemental Figure 1",
          theme = ggplot2::theme(
            plot.title.position = "plot",
            plot.title = ggplot2::element_text(hjust = 0, face = "bold", family = "Arial", size = 16),
            plot.margin = grid::unit(c(0.5, 0.5, 0.5, 0.5), "in")
          )
        )
        
    #- 6.4.3: Save Supplemental Figure 1 (8.5 x 11 with preserved aspect ratio)
      ggsave(
        "Supplemental_Figure_1_MFN_enrichment.png",
        Supplemental_Figure_1,
        width = 8.5,
        height = 11,
        units = "in",
        dpi = 300
      )      
  #+ 6.5: Figure 4 - KEGG Pathway Enrichment + Network
    #- 6.5.1: Calculate KEGG plot dimensions
      kegg_width <- length(unique(variant_enrichment_KEGG$Comparisons)) * 0.3 + 6
      kegg_height <- length(unique(variant_enrichment_KEGG$pathway_name)) * 0.3 + 2
      
    #- 6.5.2: Assemble Figure 4 (KEGG top, network bottom)
      Figure_4 <- patchwork::wrap_plots(
        variant_enrichment_plot_KEGG,   # Top half
        p,                              # Bottom half - the enrichment network plot
        nrow = 2,                       # Stacked vertically
        heights = c(1, 1)              # Equal heights for top and bottom
      ) +
        patchwork::plot_annotation(
          title = "Figure 4",
          tag_levels = "A",
          theme = ggplot2::theme(
            plot.title.position = "plot",
            plot.title = ggplot2::element_text(hjust = 0, face = "bold", family = "Arial", size = 16),
            plot.margin = grid::unit(c(0.3, 0.5, 0.3, 0.5), "in")
          )
        ) &
        ggplot2::theme(
          plot.tag.position = c(0, 0.98),
          plot.tag = ggplot2::element_text(size = 14, face = "bold", vjust = 0, hjust = 0, family = "Arial", color = "black")
        )
        
    #- 6.5.3: Save Figure 4
      ggsave(
        "Figure_4_KEGG_and_network.png",
        Figure_4,
        width = max(kegg_width, 10),    # Use wider of KEGG plot or network (10 inches)
        height = kegg_height + 8 + 1,   # KEGG height + network height (8) + title space
        units = "in",
        dpi = 300
      )
    