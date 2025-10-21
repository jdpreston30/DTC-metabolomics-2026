#* 5: Render Figures
source("R/Scripts/04_assign_plots.R")
#+ 5.1: Render Figures
print_to_png(ggdraw(xlim = c(0, 8.5), ylim = c(0, 11)) +
  draw_plot(A, x = 1.0, y = 6.60, width = 3.25, height = 3.1) +
  draw_plot(ggdraw() + draw_grob(B), x = 4.2, y = 6.62, width = 5.9 / 1.5, height = 5.4 / 1.5) +
  draw_plot(ggdraw() + draw_grob(C), x = 1, y = 3.7, width = 10.5 / 3.2, height = 10.5 / 3.2) +
  draw_plot(D.1, x = 4.75, y = 5, width = 1.5, height = 1.5) +
  draw_plot(D.2, x = 6.2, y = 5, width = 1.5, height = 1.5) +
  draw_plot(D.3, x = 4.75, y = 3.65, width = 1.5, height = 1.5) +
  draw_plot(D.4, x = 6.2, y = 3.65, width = 1.5, height = 1.5) +
  figure_labels(list(
    A = c(1.08, 9.88),
    B = c(4.40, 9.88),
    C = c(1.08, 6.43),
    D = c(4.40, 6.43)
  )), "fig1.png", width = 8.5, height = 11, dpi = 600)

#+ 5.2: Add caption to figure
source("R/Utilities/Visualization/add_caption_to_png.R")

# Define the caption text (Option 1 - simple wrapping)
figure_caption <- "Figure 1. Hierarchical clustering analysis of the top 1000 metabolic features with the greatest variance reveals to distinct clusters (A). Enrichment network plot comparing altered pathways between clusters 1 and 2 where color represents p-value, node size corresponds with enrichment factor, and connecting line thickness/node proximity indicate the number of shared pathway members (B). T-stage composition of cluster 1 versus 2 tumors (C). Annotated metabolites which showed the most significant differences between cluster 1 and 2 tumors (D)."

# Add caption to the figure
add_caption_to_png(
  png_path = "Figures/fig1.png",
  caption_text = figure_caption,
  output_path = "Figures/fig1_with_caption.png",
  caption_size = 11,
  caption_family = "Arial",
  figure_width = 8.5,
  figure_height = 11.5,  # Slightly taller to accommodate caption
  dpi = 300  # High quality
)
