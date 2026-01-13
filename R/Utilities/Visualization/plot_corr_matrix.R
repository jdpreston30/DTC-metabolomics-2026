#' Create Correlation Matrix Plot for Metabolomics Features
#'
#' This function creates a triangle correlation matrix with dots sized by p-value
#' and colored by correlation direction using corrplot.
#'
#' @param feature_table Tibble with ID column and feature columns
#' @param metadata_table Tibble with feature and display_name columns for labeling
#' @param p_threshold Numeric p-value threshold for display (default 0.1, no dot if p > threshold)#' @param show_labels Logical whether to show metabolite labels (default TRUE, set FALSE for large matrices)#' @param output_path Optional path to save as PNG (default NULL for display only)
#' @param width Width in inches for saved plot (default 10)
#' @param height Height in inches for saved plot (default 10)
#' @param dpi Resolution for saved plot (default 600)
#'
#' @return List with two elements:
#'   - summary: Overall statistics
#'     - n_metabolites: Number of metabolites
#'     - total_pairs: Total correlation pairs tested
#'     - significant_pairs: Number of significant correlations
#'     - positive_correlations: Number of significant positive correlations
#'     - negative_correlations: Number of significant negative correlations
#'     - percent_significant: Percentage of significant correlations
#'     - p_threshold: P-value threshold used
#'   - metabolite_stats: Tibble with per-metabolite statistics
#'     - feature: Feature column name
#'     - display_name: Metabolite display name
#'     - n_positive_corr: Number of significant positive correlations for this metabolite
#'     - n_negative_corr: Number of significant negative correlations for this metabolite
#'     - total_sig_corr: Total significant correlations for this metabolite
#'
#' @examples
#' \dontrun{
#' plot_corr_matrix(
#'   feature_table = TFT_cor_matrix_data,
#'   metadata_table = QC_matrix,
#'   output_path = "Outputs/Figures/corr_matrix.png"
#' )
#' }
#'
#' @export
plot_corr_matrix <- function(
    feature_table,
    metadata_table,
    p_threshold = 0.05,
    show_labels = TRUE,
    output_path = NULL,
    width = 10,
    height = 10,
    dpi = 600) {
  
  # Load required libraries
  library(dplyr)
  library(corrplot)
  
  # Remove ID column and get feature matrix
  feature_matrix <- feature_table |>
    select(-ID) |>
    as.matrix()
  
  # Calculate correlation matrix and p-values
  cor_result <- cor.mtest(feature_matrix, conf.level = 1 - p_threshold)
  cor_matrix <- cor(feature_matrix, use = "pairwise.complete.obs")
  p_matrix <- cor_result$p
  
  # Replace column names with display names
  feature_lookup <- metadata_table |>
    select(feature, display_name) |>
    deframe()
  
  display_names <- feature_lookup[colnames(cor_matrix)]
  colnames(cor_matrix) <- display_names
  rownames(cor_matrix) <- display_names
  colnames(p_matrix) <- display_names
  rownames(p_matrix) <- display_names
  
  # If output path specified, save to file
  if (!is.null(output_path)) {
    # Ensure directory exists
    output_dir <- dirname(output_path)
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    
    # Open PNG device
    png(
      filename = output_path,
      width = width,
      height = height,
      units = "in",
      res = dpi
    )
  }
  
  # Set Arial font
  par(family = "Arial")
  
  # Create the correlation plot
  corrplot(
    cor_matrix, 
    type = "lower",           # Lower triangle
    order = "hclust",         # Hierarchical clustering
    tl.col = if (show_labels) "black" else NA,  # Label color (NA hides labels)
    tl.srt = 90,              # Rotate bottom labels 90 degrees
    tl.cex = 0.8,             # Label size
    method = "circle",        # Use circles
    col = colorRampPalette(c("#113d6a", "white", "#800017"))(200),  # Blue to red
    p.mat = p_matrix,         # P-value matrix
    sig.level = p_threshold,  # Significance threshold
    insig = "blank",          # Hide non-significant
    diag = FALSE,             # Hide diagonal
    tl.pos = if (show_labels) "ld" else "n",  # Labels on left/diagonal or none
    cl.pos = "r",             # Color legend on right
    addCoef.col = NULL        # Don't add correlation coefficients
  )
  
  # Close device if saving
  if (!is.null(output_path)) {
    dev.off()
    cat("Correlation matrix saved to:", output_path, "\n")
  }
  
  # Calculate summary statistics
  total_pairs <- sum(lower.tri(cor_matrix))
  sig_pairs <- sum(p_matrix[lower.tri(p_matrix)] <= p_threshold, na.rm = TRUE)
  pos_sig <- sum(cor_matrix[lower.tri(cor_matrix)] > 0 & p_matrix[lower.tri(p_matrix)] <= p_threshold, na.rm = TRUE)
  neg_sig <- sum(cor_matrix[lower.tri(cor_matrix)] < 0 & p_matrix[lower.tri(p_matrix)] <= p_threshold, na.rm = TRUE)
  
  # Create summary list
  summary <- list(
    n_metabolites = nrow(cor_matrix),
    total_pairs = total_pairs,
    significant_pairs = sig_pairs,
    positive_correlations = pos_sig,
    negative_correlations = neg_sig,
    percent_significant = round(100 * sig_pairs / total_pairs, 1),
    p_threshold = p_threshold
  )
  
  # Calculate per-metabolite correlation statistics
  metabolite_stats <- tibble::tibble(
    feature = colnames(cor_matrix),
    display_name = display_names
  )
  
  # For each metabolite, count significant correlations
  metabolite_stats <- metabolite_stats |>
    mutate(
      n_positive_corr = purrr::map_int(feature, ~{
        idx <- which(colnames(cor_matrix) == .x)
        sum(cor_matrix[, idx] > 0 & p_matrix[, idx] <= p_threshold & 
            row.names(cor_matrix) != .x, na.rm = TRUE)
      }),
      n_negative_corr = purrr::map_int(feature, ~{
        idx <- which(colnames(cor_matrix) == .x)
        sum(cor_matrix[, idx] < 0 & p_matrix[, idx] <= p_threshold & 
            row.names(cor_matrix) != .x, na.rm = TRUE)
      }),
      total_sig_corr = n_positive_corr + n_negative_corr
    ) |>
    arrange(desc(total_sig_corr))
  
  # Print summary
  cat("\n=== Correlation Matrix Summary ===\n")
  cat("Metabolites:", summary$n_metabolites, "\n")
  cat("Total pairs tested:", summary$total_pairs, "\n")
  cat("Significant pairs (p ≤", p_threshold, "):", summary$significant_pairs, 
      "(", summary$percent_significant, "%)\n")
  cat("  - Positive correlations:", summary$positive_correlations, "\n")
  cat("  - Negative correlations:", summary$negative_correlations, "\n")
  cat("==================================\n\n")
  
  # Return both summary and per-metabolite stats
  invisible(list(
    summary = summary,
    metabolite_stats = metabolite_stats
  ))
}
