#' Simple theme for single-variable plots
#' 
#' Simplified theme perfect for PERMANOVA results, single-variable bar plots,
#' and other plots that don't need legends or complex styling. Features
#' consistent margins, panel borders, and grid lines that match PCA themes.
#' 
#' @param base_size Base font size (default: 12)
#' @param base_family Base font family (default: "Arial")
#' @return A ggplot2 theme object
#' @export
#' 
#' @examples
#' # For simple bar plots (like PERMANOVA results)
#' library(ggplot2)
#' ggplot(mtcars, aes(y = reorder(rownames(mtcars), mpg), x = mpg)) + 
#'   geom_col(fill = "black") + 
#'   theme_publication_simple() +
#'   labs(y = NULL, x = "Miles per Gallon")
#'   
#' # For vertical plots
#' ggplot(mtcars, aes(x = factor(cyl), y = mpg)) + 
#'   geom_col(fill = "black") + 
#'   theme_publication_simple()
theme_publication_simple <- function(base_size = 12, base_family = "Arial") {
  ggprism::theme_prism(base_size = base_size, base_family = base_family) +
  ggplot2::theme(
    # Plot margins - match PCA theme
    plot.margin = grid::unit(c(8, 5, 5, 5), "pt"),
    
    # Axis styling - consistent with PCA theme
    axis.line = ggplot2::element_blank(),
    axis.ticks = ggplot2::element_line(linewidth = 0.3),
    axis.ticks.length = grid::unit(0.2, "cm"),
    axis.title.y = ggplot2::element_text(size = 9, face = "bold"),
    axis.title.x = ggplot2::element_text(size = 9, face = "bold"),
    axis.text.y = ggplot2::element_text(size = 8, face = "bold"),
    axis.text.x = ggplot2::element_text(size = 8, face = "bold"),
    
    # Panel styling - match PCA border and grid
    panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.5),
    panel.grid.major.x = ggplot2::element_line(color = "gray80", linewidth = 0.3),
    panel.grid.major.y = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    panel.background = ggplot2::element_blank(),
    
    # No legend
    legend.position = "none"
  )
}
