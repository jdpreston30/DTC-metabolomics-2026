#* 6: Figure Creation
#+ 6.0: Setup
  #- 6.1.1: Setup a blank plot for spacing
    blank_plot <- ggplot2::ggplot() +
      ggplot2::theme_void() +
      ggplot2::labs(tag = NULL) + # this removes the letter completely
      ggplot2::theme(
        plot.tag = ggplot2::element_text(color = "white") # fallback if patchwork still assigns tags
      )
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
    print_to_png(Figure_1, "Figure 1")
#+ 6.2: Figure 2 - Pathway Enrichment
  #- 6.2.1: Load in KEGG and enrichment plots
    kegg_grob <- rasterGrob(as.raster(image_read("Outputs/Grob/variant_enrichment_plot_KEGG.png")), interpolate = TRUE)
    enrichment_grob <- rasterGrob(as.raster(image_read("Outputs/Grob/enrichment_network.png")), interpolate = TRUE)
  #- 6.2.2: Convert grob to ggplot for KEGG
    kegg_as_plot <- ggplot2::ggplot() +
      ggplot2::annotation_custom(
        grid::rasterGrob(as.raster(magick::image_read("Outputs/Grob/variant_enrichment_plot_KEGG.png")),
          interpolate = TRUE
        ),
        xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf
      ) +
      ggplot2::theme_void() +
      ggplot2::theme(plot.margin = grid::unit(c(0, 0, 0, 0), "in"))
  #- 6.2.4: Convert grob to ggplot for enrichment
    enrichment_as_plot <- ggplot2::ggplot() +
      ggplot2::annotation_custom(
        grid::rasterGrob(as.raster(magick::image_read("Outputs/Grob/enrichment_network.png")),
          interpolate = TRUE
        ),
        xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf
      ) +
      ggplot2::theme_void() +
      theme(plot.margin = margin(t = -0, r = 0, b = 0, l = 0))
  #- 6.2.3: Add the fixed title band via patchwork::plot_annotation
    Figure_2 <- patchwork::wrap_plots(
      kegg_as_plot + ggplot2::labs(tag = "A"),
      enrichment_as_plot + ggplot2::labs(tag = "B"),
      blank_plot + ggplot2::labs(tag = NULL),
      blank_plot + ggplot2::labs(tag = NULL),
      design = "
      AB
      CD
      "
    ) +
      patchwork::plot_annotation(
        title = "Figure 2",
        theme = ggplot2::theme(
          plot.title.position = "plot",
          plot.title = ggplot2::element_text(
            hjust = 0, face = "bold", family = "Arial", size = 16
          ),
          plot.margin = grid::unit(c(0.3, 0.5, 0.3, 0.5), "in")
        )
      ) &
      ggplot2::theme(
        plot.tag = ggplot2::element_text(
          size = 14, face = "bold", family = "Arial", color = "black"
        ),
        plot.tag.position = c(0, 0.98)
      )
  #- 6.2.4: Print
    print_to_png(Figure_2, "Figure 2", dpi = 600)
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
        plot.title = ggplot2::element_text(
          hjust = 0, face = "bold", family = "Arial", size = 16
        ),
        plot.margin = grid::unit(c(0.3, 2, 0.3, 2), "in") 
      )
    ) &
    ggplot2::theme(
      plot.tag.position = c(0, 0.98),
      plot.tag = ggplot2::element_text(size = 14, face = "bold", vjust = 0, hjust = 0, family = "Arial", color = "black")
    )
  #- 6.3.2 Save Figure 3 to PNG
    print_to_png(Figure_3, "Figure 3")
#+ 6.4: Supplemental Figure 1 - MFN Pathway Enrichment
  #- 6.4.1: Load in MFN plot as grob
    mfn_grob <- rasterGrob(as.raster(image_read("Outputs/Grob/variant_enrichment_plot_MFN.png")), interpolate = TRUE)
  #- 6.4.2: Convert grob to ggplot
    mfn_as_plot <- ggplot2::ggplot() +
      ggplot2::annotation_custom(mfn_grob, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
      ggplot2::theme_void() +
      ggplot2::theme(plot.margin = grid::unit(c(0, 0, 0, 0), "in"))
  #- 6.4.3: Add the fixed title band via patchwork::plot_annotation
    Supplemental_Figure_1 <- mfn_as_plot +
      patchwork::plot_annotation(
        title = "Supplemental Figure 1",
        theme = ggplot2::theme(
          plot.title.position = "plot",
          plot.title = ggplot2::element_text(hjust = 0, face = "bold", family = "Arial", size = 16),
          plot.margin = grid::unit(c(0.3, 0.5, 0.3, 0.5), "in") # SAME as Fig 1 & 3
        )
      )
  #- 6.4.4: Print
    print_to_png(Supplemental_Figure_1, "Supplemental Figure 1")
#+ 6.5: Build combined PDF
{
png_files <- c(
  here::here("Figures", "Figure 1.png"),
  here::here("Figures", "Figure 2.png"),
  here::here("Figures", "Figure 3.png"),
  here::here("Figures", "Supplemental Figure 1.png")
)

# Read them in
imgs <- lapply(png_files, image_read)

# Combine into a single PDF
pdf_file <- here::here("Figures", "Figures_Compiled.pdf")
image_write(image_join(imgs), path = pdf_file, format = "pdf")

message("✅ Combined PDF saved: ", pdf_file)
}
