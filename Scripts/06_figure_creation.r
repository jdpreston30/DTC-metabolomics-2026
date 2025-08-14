#* 6: Figure Creation
  #+ 6.2: Figure 1 - PCA Analysis (4 columns x 3 rows)
    #- 6.2.1: Create PCA panels with small points (using existing data and colors)
      pca_panel_A <- make_PCA(
        data = TFT_metaboanalyst_log2_path %>% select(-c(T_computed_bin:LVI)),
        ellipse_colors = variant_colors,
        point_size = 0.5
      )$plot
      pca_panel_B <- make_PCA(
        data = TFT_metaboanalyst_log2_path %>% 
          select(-c(LD:Variant)) %>%
          rename(Variant = T_computed_bin) %>%
          filter(!is.na(Variant)),
        ellipse_colors = T_stage_colors,
        point_size = 0.5
      )$plot
      pca_panel_C <- make_PCA(
        data = TFT_metaboanalyst_log2_path %>% 
          select(-c(T_computed_bin, LD, Variant)) %>%
          mutate(LVI = ifelse(is.na(LVI), "NA", as.character(LVI))) %>% 
          rename(Variant = LVI),
        ellipse_colors = LVI_colors,
        point_size = 0.5
      )$plot
      pca_panel_D <- permanova_table_plot  # Replace placeholder with PERMANOVA table
    #- 6.2.2: Apply consistent theming for square PCAs
      pca_theme <- ggplot2::theme(
        plot.margin = grid::unit(c(2, 2, 2, 2), "pt"),
        aspect.ratio = 1,  # Perfect 1:1 aspect ratio for squares
        # Override the huge text sizes from original PCA plots
        axis.title = ggplot2::element_text(size = 8, face = "bold"),
        axis.text = ggplot2::element_text(size = 7, face = "bold", color = "black"),
        panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.5),
        panel.grid.major = ggplot2::element_line(color = "gray80", linewidth = 0.3),
        # Add legend above each plot, horizontal layout
        legend.position = "top",
        legend.direction = "horizontal",
        legend.box = "horizontal",
        legend.margin = ggplot2::margin(0, 0, 5, 0),
        legend.text = ggplot2::element_text(size = 7),
        legend.title = ggplot2::element_blank(),  # Remove "Class" label
        legend.key.size = grid::unit(0.3, "cm")
      )
      
      # Apply theme only - no extra geom_point layers needed
      pca_A_square <- pca_panel_A + pca_theme
      pca_B_square <- pca_panel_B + pca_theme  
      pca_C_square <- pca_panel_C + pca_theme
    #- 6.2.3: Create heatmap plot for bottom row
      # Convert pheatmap object to grob for patchwork compatibility
      heatmap_bottom <- patchwork::wrap_elements(variance_500$heatmap_plot$gtable)
    #- 6.2.4: Assemble Figure 1 (3 columns top row + full width bottom row)
      Figure_1 <- patchwork::wrap_plots(
        # Row 1: 3 PCA plots (A, B, C)
        pca_A_square, pca_B_square, pca_C_square,
        # Row 2: 1 heatmap spanning full width (D)
        heatmap_bottom,
        design = "ABC\nDDD",  # Top row: 3 equal columns, Bottom row: spans all 3 columns
        heights = c(1, 3)  # PCA row: 1/4 of page, Heatmap: 3/4 of page
      ) +
        patchwork::plot_annotation(
          title = "Figure 1: Principal Component Analysis",
          tag_levels = "A",
          theme = ggplot2::theme(
            plot.title.position = "plot",
            plot.title = ggplot2::element_text(hjust = 0, face = "bold", family = "Arial", size = 15),
            plot.margin = grid::unit(c(0.3, 0.5, 0.3, 0.5), "in"),
            plot.tag.position = c(0, 1),
            plot.tag = ggplot2::element_text(size = 14, face = "bold", vjust = 0, hjust = 0)
          )
        )
    #- 6.2.5: Preview Figure 1
      smart_preview(Figure_1)

  #+ 6.3: Figure 2 - PERMANOVA Results Table (2x2 Layout)
    #- 6.3.1: Create placeholder plots for positions B, C, D
      placeholder_plot <- ggplot2::ggplot() + 
        ggplot2::theme_void() +
        ggplot2::labs(title = "Placeholder") +
        ggplot2::theme(
          plot.title = ggplot2::element_text(hjust = 0.5, size = 12, color = "gray50")
        )
    #- 6.3.2: Create Figure 2 with 2x2 layout (table in top-left)
      Figure_2 <- patchwork::wrap_plots(
        # Row 1: PERMANOVA table (A) + placeholder (B)
        permanova_table_plot, placeholder_plot,
        # Row 2: placeholder (C) + placeholder (D)  
        placeholder_plot, placeholder_plot,
        ncol = 2, 
        nrow = 2
      ) +
        patchwork::plot_annotation(
          title = "Figure 2: PERMANOVA Analysis Results",
          tag_levels = "A",
          theme = ggplot2::theme(
            plot.title.position = "plot",
            plot.title = ggplot2::element_text(hjust = 0, face = "bold", family = "Arial", size = 15),
            plot.margin = grid::unit(c(0.3, 0.5, 0.3, 0.5), "in"),
            plot.tag.position = c(0, 1),
            plot.tag = ggplot2::element_text(size = 14, face = "bold", vjust = 0, hjust = 0)
          )
        )
    #- 6.3.3: Preview Figure 2
      smart_preview(Figure_2)

  #+ 6.4: Export High-Resolution PNGs
    #- 6.4.1: Create directory if needed
      dir.create("Outputs", showWarnings = FALSE)
    #- 6.4.2: Export Figure 1 as high-res PNG (8.5 x 11 @ 600 DPI)
      ggplot2::ggsave(
        "Outputs/Figure_1.png",
        Figure_1,
        width = 8.5, height = 11, units = "in",
        dpi = 600, bg = "white"
      )
      cat("Figure 1 exported to: Outputs/Figure_1.png\n")
    #- 6.4.3: Export Figure 2 as high-res PNG (8.5 x 11 @ 600 DPI)  
      ggplot2::ggsave(
        "Outputs/Figure_2.png",
        Figure_2,
        width = 8.5, height = 11, units = "in",
        dpi = 600, bg = "white"
      )
      cat("Figure 2 exported to: Outputs/Figure_2.png\n")

  #+ 6.5: Figure 3 - Targeted Metabolomics 4-up
    #- 6.2.1: Define individual panels
      # For now using same plot for all panels - replace as needed
      fig3_panel_A <- fig3_A
      fig3_panel_B <- fig3_A  # Replace with actual plot B
      fig3_panel_C <- fig3_A  # Replace with actual plot C  
      fig3_panel_D <- fig3_A  # Replace with actual plot D
    #- 6.2.2: Apply consistent theming to panels
      panel_theme <- ggplot2::theme(
        plot.margin = grid::unit(c(4, 4, 4, 4), "pt"),
        aspect.ratio = 2.5/3  # height/width ratio
      )
      
      fig3_A_small <- fig3_panel_A + panel_theme
    #- 6.2.3: Assemble Figure 3 (Large top plot + two bottom plots)
      Figure_3 <- patchwork::wrap_plots(
        # Top row: single plot spanning full width
        fig3_A_small,
        # Bottom row: two plots side by side  
        fig3_C_small, fig3_D_small,
        # Layout specification
        design = "AA\nBC",  # A spans both columns on top, B and C on bottom
        heights = c(1, 1)  # Equal height for top and bottom rows
      ) +
        patchwork::plot_annotation(
          title = "Figure 3: Targeted Metabolomics Analysis",
          tag_levels = "A",
          theme = ggplot2::theme(
            plot.title.position = "plot",
            plot.title = ggplot2::element_text(hjust = 0, face = "bold", family = "Arial", size = 15),
            plot.margin = grid::unit(c(0.3, 0.5, 0.3, 0.5), "in"),
            plot.tag.position = c(0, 1),
            plot.tag = ggplot2::element_text(size = 14, face = "bold", vjust = 0, hjust = 0)
          )
        )
    #- 6.2.4: Preview Figure 3
      smart_preview(Figure_3)



  #+ 6.5: Export functions (uncomment when ready)
    #- 6.5.1: High-resolution exports
      # ggsave("Outputs/Figure_3.pdf", Figure_3, width = 8.5, height = 11, dpi = 300)
      # ggsave("Outputs/Figure_3.png", Figure_3, width = 8.5, height = 11, dpi = 300)
