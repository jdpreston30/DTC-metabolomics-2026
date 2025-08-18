#' Dot plot theme with floating bars for publication
#' 
#' Theme specifically designed for dot plots with floating mean bars,
#' includes legend positioning for grouped data. Perfect for showing
#' individual data points with group means.
#' 
#' @param base_size Base font size (default: 12)
#' @param base_family Base font family (default: "Arial")
#' @param legend_position Position of legend (default: c(0.05, 0.8))
#' @return A ggplot2 theme object
#' @export
#' 
#' @examples
#' # For grouped dot plots with legend
#' library(ggplot2)
#' ggplot(mtcars, aes(x = factor(cyl), y = mpg, fill = factor(gear))) + 
#'   geom_col(position = "dodge", alpha = 0.7) + 
#'   geom_jitter(width = 0.2) +
#'   theme_publication_dotplot()
#'   
#' # Custom legend position
#' ggplot(mtcars, aes(x = factor(cyl), y = mpg, fill = factor(gear))) + 
#'   geom_col(position = "dodge", alpha = 0.7) + 
#'   theme_publication_dotplot(legend_position = c(0.8, 0.8))
theme_publication_dotplot <- function(base_size = 12, base_family = "Arial", 
                                     legend_position = "top", legend_direction = "horizontal") {
  ggprism::theme_prism(base_size = base_size, base_family = base_family) +
  ggplot2::theme(
    # Plot margins
    plot.margin = grid::unit(c(8, 5, 5, 5), "pt"),
    # Axis styling
    axis.line = ggplot2::element_line(linewidth = 0.5, lineend = "square"),
    axis.ticks = ggplot2::element_line(linewidth = 0.3),
    axis.ticks.length = grid::unit(0.2, "cm"),
    axis.title.y = ggplot2::element_text(size = 9, face = "bold"),
    axis.title.x = ggplot2::element_text(size = 9, face = "bold"),
    axis.text.y = ggplot2::element_text(size = 8, face = "bold"),
    axis.text.x = ggplot2::element_text(size = 8, face = "bold", angle = 45, hjust = 1, vjust = 0.97),
    # Panel styling
    panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.5),
    panel.grid.major.x = ggplot2::element_line(color = "gray80", linewidth = 0.3),
    panel.grid.major.y = ggplot2::element_blank(),
    panel.grid.minor.x = ggplot2::element_blank(),
    panel.grid.minor.y = ggplot2::element_blank(),
    panel.background = ggplot2::element_blank(),
    panel.ontop = FALSE,
    # Legend styling for grouped plots
    legend.title = ggplot2::element_blank(),
    legend.position = legend_position,
    legend.direction = legend_direction,
    legend.justification = "center",
    legend.key.size = grid::unit(0.4, "cm"),
    legend.key.width = grid::unit(0.4, "cm"),
    legend.key.height = grid::unit(0.4, "cm"),
    legend.text = ggplot2::element_text(size = 8, face = "bold")
  )
}
