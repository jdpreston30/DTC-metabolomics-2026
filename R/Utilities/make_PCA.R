#' Create and save PCA plot with ellipses
#'
#' @param data Data frame with Patient_ID, Variant, and feature columns
#' @param output_filename Filename for the saved SVG plot (without path)
#' @param plot_title Optional title for the plot (default: "")
#' @param ellipse_colors Named vector of colors for each variant/group
#' @return List containing the plot, PCA object, scores, scores_df, and explained variance
#' @export
make_PCA <- function(data, output_filename, plot_title = "", 
                        ellipse_colors = c("PTC" = "#DF8D0A", "FV-PTC" = "#23744E", "FTC" = "#194992")) {
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
          size = 3, shape = 21, stroke = 0.8,
          fill = ellipse_colors[as.character(scores_df$Class)]
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
      
      #_Save plot
      ggplot2::ggsave(
        filename = paste0("Outputs/", output_filename),
        plot = pca_plot,
        device = "svg",
        width = 8,
        height = 8,
        units = "in",
        dpi = 600
      )
      
      #_Print plot and return results
      print(pca_plot)
      cat("PCA plot saved to:", paste0("Outputs/", output_filename), "\n")
      cat("PC1 explains", explained[1], "% of variance\n")
      cat("PC2 explains", explained[2], "% of variance\n")
      
      #_Return useful objects for further analysis
      return(list(
        plot = pca_plot,
        pca = pca,
        scores = scores,
        scores_df = scores_df,
        explained = explained
      ))
}
