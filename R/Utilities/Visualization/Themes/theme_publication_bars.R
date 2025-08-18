#' Standard bar plot theme for publication
#' 
#' A clean, professional theme based on ggprism for bar plots and similar
#' visualizations. Uses Arial font and consistent styling with angled x-axis labels.
#' 
#' @param base_size Base font size (default: 12)
#' @param base_family Base font family (default: "Arial")
#' @return A ggplot2 theme object
#' @export
#' 
#' @examples
#' # Basic usage
#' library(ggplot2)
#' ggplot(mtcars, aes(x = factor(cyl), y = mpg)) + 
#'   geom_col() + 
#'   theme_publication_bars()
#'   
#' # With custom fill colors
#' ggplot(mtcars, aes(x = factor(cyl), y = mpg)) + 
#'   geom_col(fill = "black") + 
#'   theme_publication_bars()
theme_publication_bars <- function(base_size = 12, base_family = "Arial") {
  ggprism::theme_prism(base_size = base_size, base_family = base_family) +
  ggplot2::theme(
    # Axis styling
    axis.line = ggplot2::element_line(linewidth = 1.5, lineend = "square"),
    axis.ticks = ggplot2::element_line(linewidth = 1),
    axis.ticks.length = grid::unit(0.2, "cm"),
    axis.title.y = ggplot2::element_text(size = base_size, face = "bold"),
    axis.text.y = ggplot2::element_text(size = base_size - 1, face = "bold"),
    axis.text.x = ggplot2::element_text(size = base_size - 2, face = "bold", 
                                       angle = 45, hjust = 1, vjust = 0.97),
    
    # Panel styling
    panel.grid = ggplot2::element_blank(),
    
    # Legend styling
    legend.title = ggplot2::element_blank(),
    legend.position = "none"
  )
}
