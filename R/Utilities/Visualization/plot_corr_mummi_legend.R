#' Create Correlation Mummichog-Style Legend
#'
#' Creates a standalone legend for -log10(p-value) correlation enrichment results
#' with a gradient from dark blue (0.6) to dark red (2+). The asterisk at 1.30103
#' (p=0.05) marks the color transition midpoint.
#'
#' @return A grob object containing only the legend (no plot)
#'
#' @examples
#' \dontrun{
#' corr_mummi_legend <- plot_corr_mummi_legend()
#' }
#'
#' @export
plot_corr_mummi_legend <- function() {
  
  # Load required libraries
  library(ggplot2)
  library(cowplot)
  
  # Create dummy data for the legend
  legend_data <- data.frame(
    x = 1,
    y = seq(0.6, 2, length.out = 100),
    p_value = seq(0.6, 2, length.out = 100)
  )
  
  # Breaks: 0.6, 1, asterisk (1.30103), 1.5, 2+
  asterisk_log <- -log10(0.05)
  all_breaks <- c(0.6, 1, asterisk_log, 1.5, 2)
  all_labels <- c("0.6", "1", "✱", "1.5", "2+")
  
  # Create the plot
  p <- ggplot(legend_data, aes(x = x, y = y, fill = p_value)) +
    geom_tile() +
    scale_fill_gradientn(
      colors = c("#113d6a", "#D8919A", "#800017"),
      values = c(0, (asterisk_log - 0.6) / (2 - 0.6), 1),
      limits = c(0.6, 2),
      breaks = all_breaks,
      labels = all_labels,
      name = "-log10(p-value)"
    ) +
    guides(
      fill = guide_colorbar(
        barheight = unit((8.22 - 0.11853333333333334) * 0.9, "cm"),
        barwidth = unit(0.27, "cm"),
        ticks.colour = NA,
        ticks.linewidth = 0,
        direction = "vertical",
        label.position = "right",
        title.position = "left",
        frame.colour = "black",
        frame.linewidth = 0.5,
        draw.ulim = TRUE,
        draw.llim = TRUE
      )
    ) +
    theme_void(base_family = "Arial") +
    theme(
      text = element_text(family = "Arial"),
      legend.title = element_text(size = 11, face = "bold", family = "Arial", angle = 90, vjust = 0.5, hjust = 0.5),
      legend.text = element_text(size = 11, family = "Arial"),
      legend.position = "right"
    )
  
  # Extract just the legend as a grob
  legend <- get_legend(p)
  
  # Return only the legend
  return(legend)
}
