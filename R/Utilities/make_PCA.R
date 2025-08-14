#' Create PCA plot with ellipses
#'
#' @param data Data frame with Patient_ID, Variant, and feature columns
#' @param plot_title Optional title for the plot (default: "")
#' @param ellipse_colors Named vector of colors for each variant/group
#' @param point_size Size of the points (default: 3 for standalone, 0.5 for multi-panel)
#' @return List containing the plot, PCA object, scores, scores_df, and explained variance
#' @export
make_PCA <- function(data, plot_title = "", 
                        ellipse_colors = c("PTC" = "#DF8D0A", "FV-PTC" = "#23744E", "FTC" = "#194992"),
                        point_size = 3) {
      #_Data preparation
      df <- as.data.frame(data)
      cls_col <- if ("Variant" %in% names(df)) "Variant" else names(df)[2]
      X <- df[, -c(1, 2), drop = FALSE]
      #_Coerce to numeric safely
      X[] <- lapply(X, function(v) suppressWarnings(as.numeric(v)))
      #_Handle NAs with median imputation
      if (anyNA(X)) {
        X[] <- lapply(X, function(v) {
          v[is.na(v)] <- stats::median(v, na.rm = TRUE)
          v
        })
      }
      Y <- factor(df[[cls_col]])
      
      #_Perform PCA
      pca <- stats::prcomp(X, center = TRUE, scale. = TRUE)
      scores <- pca$x[, 1:2, drop = FALSE]
      explained <- round((pca$sdev^2 / sum(pca$sdev^2))[1:2] * 100)
      
      #_Prepare plot data
      scores_df <- data.frame(
        Comp1 = scores[, 1],
        Comp2 = scores[, 2],
        Class = Y
      )
      
      #_Create PCA plot
      pca_plot <- ggplot2::ggplot(scores_df, ggplot2::aes(x = Comp1, y = Comp2, color = Class)) +
        ggplot2::geom_point(
          size = point_size, shape = 16  # Use parameter for size, solid circles
        ) +
        ggplot2::stat_ellipse(
          geom = "polygon", ggplot2::aes(fill = Class),
          alpha = 0.3, color = NA
        ) +
        ggplot2::scale_color_manual(values = ellipse_colors, drop = FALSE) +
        ggplot2::scale_fill_manual(values = ellipse_colors, drop = FALSE) +
        ggplot2::theme_minimal(base_family = "Arial") +
        ggplot2::labs(
          x = paste0("PC1 (", explained[1], "%)"),
          y = paste0("PC2 (", explained[2], "%)")
        ) +
        ggplot2::theme(
          axis.title = ggplot2::element_text(size = 25, face = "bold"),
          axis.text = ggplot2::element_text(size = 22, face = "bold", color = "black"),
          legend.position = "none",
          panel.grid.major = ggplot2::element_line(color = "gray80", linewidth = 0.8, linetype = "solid"),
          panel.grid.minor = ggplot2::element_blank(),
          panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 3.2),
          panel.background = ggplot2::element_blank()
        )
      
      #_Return useful objects for further analysis
      return(list(
        plot = pca_plot,
        pca = pca,
        scores = scores,
        scores_df = scores_df,
        explained = explained
      ))
}
