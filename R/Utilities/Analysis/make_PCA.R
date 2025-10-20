#' Create PCA or PLS-DA plot with ellipses
#'
#' @param data Data frame with grouping variable and feature columns
#' @param group_var Column name to use for grouping (default: "Variant")
#' @param method Method to use: "PCA" or "PLSDA" (default: "PCA")
#' @param ncomp Number of components for PLS-DA (default: 2)
#' @param plot_title Optional title for the plot (default: "")
#' @param ellipse_colors Named vector of colors for each variant/group
#' @param point_size Size of the points (default: 3 for standalone, 0.5 for multi-panel)
#' @param comp_x Which principal component to plot on x-axis (default: 1)
#' @param comp_y Which principal component to plot on y-axis (default: 2)
#' @return List containing the plot, model object, scores, scores_df, and explained variance
#' @export
make_PCA <- function(data, group_var = "Variant", method = "PCA", ncomp = 2, plot_title = "", 
                        ellipse_colors = c("PTC" = "#DF8D0A", "FV-PTC" = "#23744E", "FTC" = "#194992"),
                        point_size = 3, comp_x = 1, comp_y = 2) {
      #_Data preparation
      df <- as.data.frame(data)
      
      # Validate group_var exists
      if (!group_var %in% names(df)) {
        stop(paste("Group variable", group_var, "not found in data"))
      }
      
      # Select ONLY metabolomics feature columns (C18 and HILIC)
      feature_pattern <- "^(C18|HILIC)"
      feature_cols <- names(df)[grepl(feature_pattern, names(df))]
      
      # Check if we have any feature columns
      if (length(feature_cols) == 0) {
        stop("No metabolomics feature columns (C18/HILIC) found in data")
      }
      
      # Create clean dataset with only features + grouping variable
      clean_data <- df[, c(group_var, feature_cols), drop = FALSE]
      X <- clean_data[, feature_cols, drop = FALSE]
      Y <- factor(clean_data[[group_var]])
      #_Coerce to numeric safely
      X[] <- lapply(X, function(v) suppressWarnings(as.numeric(v)))
      #_Handle NAs with median imputation
      if (anyNA(X)) {
        X[] <- lapply(X, function(v) {
          v[is.na(v)] <- stats::median(v, na.rm = TRUE)
          v
        })
      }
      
      #_Remove zero-variance columns
      zero_var_cols <- sapply(X, function(col) {
        var_val <- var(col, na.rm = TRUE)
        is.na(var_val) || var_val == 0
      })
      if (any(zero_var_cols)) {
        removed_cols <- names(zero_var_cols)[zero_var_cols]
        cat("Removing", length(removed_cols), "zero-variance columns:\n")
        cat(paste(removed_cols, collapse = ", "), "\n")
        X <- X[, !zero_var_cols, drop = FALSE]
      }
      
      # Check if we have any columns left
      if (ncol(X) == 0) {
        stop("No valid feature columns remaining after filtering")
      }
      
      Y <- factor(df[[group_var]])
      
      #_Perform analysis based on method
      method <- match.arg(method, c("PCA", "PLSDA"))
      
      if (method == "PCA") {
        model <- stats::prcomp(X, center = TRUE, scale. = TRUE)
        max_comp <- min(ncol(X), nrow(X) - 1)
        if (comp_x > max_comp || comp_y > max_comp) {
          stop(paste("Requested components exceed available components. Max components:", max_comp))
        }
        scores <- model$x[, c(comp_x, comp_y), drop = FALSE]
        explained <- round((model$sdev^2 / sum(model$sdev^2))[c(comp_x, comp_y)] * 100)
        comp_label <- "PC"
      } else {
        # PLS-DA using mixOmics
        if (!requireNamespace("mixOmics", quietly = TRUE)) {
          stop("Package 'mixOmics' is required for PLS-DA. Please install it.")
        }
        model <- mixOmics::plsda(X, Y, ncomp = ncomp)
        max_comp <- min(ncomp, ncol(X), nrow(X) - 1)
        if (comp_x > max_comp || comp_y > max_comp) {
          stop(paste("Requested components exceed available components. Max components:", max_comp))
        }
        scores <- model$variates$X[, c(comp_x, comp_y), drop = FALSE]
        explained <- round(model$prop_expl_var$X[c(comp_x, comp_y)] * 100)
        comp_label <- "LV"
      }
      
      #_Prepare plot data
      scores_df <- data.frame(
        Comp1 = scores[, 1],
        Comp2 = scores[, 2],
        Class = Y
      )
      
      # Identify NA values and create separate datasets
      na_mask <- is.na(scores_df$Class)
      scores_df_complete <- scores_df[!na_mask, , drop = FALSE]
      scores_df_na <- scores_df[na_mask, , drop = FALSE]
      
      #_Create PCA plot
      pca_plot <- ggplot2::ggplot() +
        # Plot complete cases with colors and ellipses
        {if(nrow(scores_df_complete) > 0) {
          list(
            ggplot2::geom_point(
              data = scores_df_complete,
              ggplot2::aes(x = Comp1, y = Comp2, color = Class, fill = Class),
              size = point_size, shape = 16
            ),
            ggplot2::stat_ellipse(
              data = scores_df_complete,
              ggplot2::aes(x = Comp1, y = Comp2, fill = Class),
              geom = "polygon", alpha = 0.3, color = NA
            )
          )
        }} +
        # Plot NA values as open circles without color
        {if(nrow(scores_df_na) > 0) {
          ggplot2::geom_point(
            data = scores_df_na,
            ggplot2::aes(x = Comp1, y = Comp2),
            size = point_size, shape = 1, color = "black", fill = NA
          )
        }} +
        ggplot2::scale_color_manual(values = ellipse_colors, drop = TRUE, na.translate = FALSE) +
        ggplot2::scale_fill_manual(values = ellipse_colors, drop = TRUE, na.translate = FALSE) +
        ggplot2::theme_minimal(base_family = "Arial") +
        ggplot2::labs(
          x = paste0(comp_label, comp_x, " (", explained[1], "%)"),
          y = paste0(comp_label, comp_y, " (", explained[2], "%)")
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
        method = method,
        plot = pca_plot,
        model = model,
        scores = scores,
        scores_df = scores_df,
        explained = explained,
        comp_label = comp_label
      ))
}
