#* 6: Figure Creation
  #+ 6.2: Figure 1 - Clustering Analysis Layout
    #- 6.2.4: Assemble Figure 1 (Heatmap top 50%, vertical PERMANOVA + 4 PCAs bottom 50%)
      Figure_1 <- patchwork::wrap_plots(
        # Top 50%: Heatmap spanning full width
        heatmap_1A,
        # Bottom 50%: Vertical PERMANOVA (left) + 4 PCAs (2x2 grid, right)
        permanova_1B, pca_1C, pca_1D, pca_1F,
        pca_1E, ,
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
    #- 6.2.6: Save Figure 1 to PDF and PNG
      print_to_png(Figure_1, "Figure_1_preview")
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
            plot.tag = ggplot2::element_text(size = 16, face = "bold", vjust = 0, hjust = 0, family = "Arial")
          )
        )
    #- 6.2.4: Preview Figure 3
      smart_preview(Figure_3)



  #+ 6.5: Export functions (uncomment when ready)
    #- 6.5.1: High-resolution exports
      # ggsave("Outputs/Figure_3.pdf", Figure_3, width = 8.5, height = 11, dpi = 300)
      # ggsave("Outputs/Figure_3.png", Figure_3, width = 8.5, height = 11, dpi = 300)
