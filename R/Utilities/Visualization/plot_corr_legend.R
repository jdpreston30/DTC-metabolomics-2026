#' Create Correlation Legend
#'
#' Creates a standalone legend for correlation values with a gradient from -1 (blue)
#' to 0 (white) to 1 (red), styled to match other publication figures.
#'
#' @return A grob object containing only the legend (no plot)
#'
#' @examples
#' \dontrun{
#' corr_legend <- plot_corr_legend()
#' }
#'
#' @export
plot_corr_legend <- function() {
  
  # Load required libraries
  library(ggplot2)
  library(cowplot)
  
  # Create dummy data for the legend
  legend_data <- data.frame(
    x = 1,
    y = seq(-1, 1, length.out = 100),
    correlation = seq(-1, 1, length.out = 100)
  )
  
  # Create the plot
  p <- ggplot(legend_data, aes(x = x, y = y, fill = correlation)) +
    geom_tile() +
    scale_fill_gradientn(
      colors = c("#113d6a", "white", "#800017"),
      limits = c(-1, 1),
      breaks = seq(-1, 1, by = 0.2),
      labels = c("-1", "-0.8", "-0.6", "-0.4", "-0.2", "0", "0.2", "0.4", "0.6", "0.8", "1"),
      name = "Pearson\nCorrelation"
    ) +
    guides(
      fill = guide_colorbar(
        barheight = unit(8.22-0.11853333333333334, "cm"),
        barwidth = unit(0.27, "cm"),
        ticks.colour = NA,
        ticks.linewidth = 0,
        direction = "vertical",
        label.position = "right",
        frame.colour = "black",
        frame.linewidth = 0.5,
        draw.ulim = TRUE,
        draw.llim = TRUE
      )
    ) +
    theme_void(base_family = "Arial") +
    theme(
      text = element_text(family = "Arial"),
      legend.title = element_text(size = 11, face = "bold", family = "Arial", vjust = 0.3, margin = margin(b = 10)),
      legend.text = element_text(size = 11, family = "Arial"),
      legend.position = "right"
    )
  
  # Extract just the legend as a grob
  legend <- get_legend(p)
  
  # Return only the legend
  return(legend)
}
