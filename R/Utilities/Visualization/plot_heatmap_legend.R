#' Create Heatmap Legend (Horizontal)
#'
#' Creates a horizontal legend for heatmap z-scores with colors matching
#' the heatmap's color scale. Styled to match other publication figures.
#' All sizing is controlled by draw_plot() width and height parameters.
#'
#' @param legend_params List containing color scale parameters from make_heatmap:
#'   - colors: Vector of colors from the heatmap
#'   - limits: c(min, max) values for the scale
#'   - breaks: Numeric vector of break points
#'   - labels: Character vector of labels for breaks
#'   - title: Legend title (default "Z-Score")
#'
#' @return A grob object containing only the legend (no plot)
#'
#' @examples
#' \dontrun{
#' heatmap_result <- make_heatmap(data, ...)
#' heatmap_legend <- plot_heatmap_legend(heatmap_result$legend_params)
#' # Control size in draw_plot: draw_plot(heatmap_legend, width = 4, height = 1)
#' }
#'
#' @export
plot_heatmap_legend <- function(legend_params) {
  
  # Load required libraries
  library(ggplot2)
  library(cowplot)
  
  # Extract parameters
  colors <- legend_params$colors
  limits <- legend_params$limits
  breaks <- legend_params$breaks
  labels <- legend_params$labels
  title <- legend_params$title
  fontsize <- legend_params$fontsize
  
  # Create dummy data for the legend
  legend_data <- data.frame(
    x = 1,
    y = seq(limits[1], limits[2], length.out = length(colors)),
    zscore = seq(limits[1], limits[2], length.out = length(colors))
  )
  
  # Create the plot - bar dimensions match draw_plot() width/height exactly
  p <- ggplot(legend_data, aes(x = x, y = y, fill = zscore)) +
    geom_tile() +
    scale_fill_gradientn(
      colors = colors,
      limits = limits,
      breaks = breaks,
      labels = labels,
      name = NULL
    ) +
    guides(
      fill = guide_colorbar(
        barheight = unit(40/300, "inches"),
        barwidth = unit(602/300, "inches"),
        ticks.colour = NA,
        ticks.linewidth = 0,
        direction = "horizontal",
        label.position = "bottom",
        label.vjust = 1.5,  # Reduce whitespace
        title.position = "top",
        frame.colour = NA,
        frame.linewidth = 0,
        draw.ulim = FALSE,
        draw.llim = FALSE,
        nbin = 300
      )
    ) +
    theme_void(base_family = "Arial") +
    theme(
      text = element_text(family = "Arial"),
      legend.title = element_text(size = 6.55, face = "bold", family = "Arial", angle = 0, vjust = 0.5, hjust = 0.5),
      legend.text = element_text(size = 6.55, family = "Arial"),
      legend.position = "bottom"
    )
  
  # Extract just the legend as a grob
  legend <- get_legend(p)
  
  # Return only the legend
  return(legend)
}
