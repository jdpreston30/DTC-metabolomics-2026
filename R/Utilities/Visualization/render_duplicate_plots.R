render_duplicate_plots <- function(duplicate_plots, output_path, feature_names) {
  #' Render Duplicate Feature Plots to PDF
  #' 
  #' Creates a PDF with duplicate metabolite features arranged in pairs for comparison.
  #' Plots are arranged 3 pairs per page (3 rows x 2 columns) on 8.5 x 11 inch pages.
  #' 
  #' @param duplicate_plots Named list of ggplot objects from plot_stage_targeted()
  #' @param output_path Path to save the PDF file
  #' @param feature_names Vector of feature column names in same order as duplicate_plots
  #' 
  #' @return Invisibly returns the output path
  #' @export
  
  # Get plot names and group them
  plot_names <- names(duplicate_plots)
  unique_names <- unique(plot_names)
  
  # Prepare all pairs with feature names as subtitles
  all_pairs <- list()
  for (name in unique_names) {
    pair_indices <- which(plot_names == name)
    pair_plots <- duplicate_plots[pair_indices]
    pair_features <- feature_names[pair_indices]
    
    # Add feature names as subtitles
    pair_plots_labeled <- purrr::map2(pair_plots, pair_features, ~{
      .x + ggplot2::labs(subtitle = .y)
    })
    
    all_pairs <- c(all_pairs, pair_plots_labeled)
  }
  
  # Calculate number of pages needed (6 plots per page = 3 pairs per page)
  n_plots <- length(all_pairs)
  plots_per_page <- 6
  n_pages <- ceiling(n_plots / plots_per_page)
  
  # Create PDF with 8.5 x 11 dimensions
  grDevices::cairo_pdf(output_path, width = 8.5, height = 11)
  
  for (page in 1:n_pages) {
    # Get plots for this page
    start_idx <- (page - 1) * plots_per_page + 1
    end_idx <- min(page * plots_per_page, n_plots)
    page_plots <- all_pairs[start_idx:end_idx]
    
    # Create grid (3 rows x 2 columns)
    combined <- cowplot::plot_grid(
      plotlist = page_plots, 
      ncol = 2, 
      nrow = 3
    )
    print(combined)
  }
  
  grDevices::dev.off()
  
  message("Duplicate plots rendered to: ", output_path)
  invisible(output_path)
}
