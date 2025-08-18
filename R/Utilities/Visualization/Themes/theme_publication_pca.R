#' PCA theme for publication figures
#' 
#' Specialized theme for PCA plots with square aspect ratio, minimal styling,
#' and horizontal legend positioning. Designed for multi-panel PCA figures
#' with consistent margins and close axis labels.
#' 
#' @param base_size Base font size (default: 8)
#' @param base_family Base font family (default: "Arial")
#' @return A ggplot2 theme object
#' @export
#' 
#' @examples
#' # For PCA plots with groups
#' library(ggplot2)
#' # Simulate PCA data
#' pca_data <- data.frame(
#'   PC1 = rnorm(100),
#'   PC2 = rnorm(100),
#'   group = rep(c("A", "B", "C"), length.out = 100)
#' )
#' 
#' ggplot(pca_data, aes(PC1, PC2, color = group)) + 
#'   geom_point() + 
#'   theme_publication_pca() +
#'   labs(x = "PC1 (45.2%)", y = "PC2 (23.1%)")
#'   
#' # With ellipses
#' ggplot(pca_data, aes(PC1, PC2, color = group)) + 
#'   geom_point() +
#'   stat_ellipse() +
#'   theme_publication_pca()
theme_publication_pca <- function(base_size = 8, base_family = "Arial") {
  ggplot2::theme(
    # Plot margins and aspect ratio
    plot.margin = grid::unit(c(2, 8, 8, 8), "pt"),
    aspect.ratio = 1,
    
    # Axis styling - closer labels
    axis.title.x = ggplot2::element_text(size = 9, face = "bold", 
                                        margin = ggplot2::margin(t = 2, r = 0, b = 0, l = 0)),
    axis.title.y = ggplot2::element_text(size = 9, face = "bold", 
                                        margin = ggplot2::margin(t = 0, r = -2, b = 0, l = 0)),
    axis.text = ggplot2::element_text(size = 8, face = "bold", color = "black"),
    
    # Panel styling
    panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.5),
    panel.grid.major = ggplot2::element_line(color = "gray80", linewidth = 0.3),
    panel.background = ggplot2::element_blank(),
    
    # Legend styling - horizontal at top
    legend.position = "top",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.margin = ggplot2::margin(t = 0, r = 0, b = -10, l = 0),
    legend.text = ggplot2::element_text(size = 8, face = "bold"),
    legend.title = ggplot2::element_blank(),
    legend.key.size = grid::unit(0.3, "cm")
  )
}
