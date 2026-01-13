#' Create Multi-Column Pathway Enrichment Dot Plot from Enrichment Data
#'
#' This function creates a dot plot visualization for multiple pathway enrichment sets
#' arranged in columns. Based on plot_mummichog_enrichment but accepts pre-processed
#' enrichment data with a column indicator for faceting.
#'
#' @param enrichment_data Data frame with columns: pathway_name, Comparisons, enrichment_factor, p_fisher
#' @param p_display_threshold Numeric p-value threshold for filtering pathways (default 0.05)
#' @param p_color Numeric vector of p-values defining color gradient (default c(0.1, 0.05, 0.01)). First value maps to darkest blue, last to darkest red.
#' @param detect_thresh Numeric minimum number of comparisons (metabolites) a pathway must be detected in to be included (default 1)
#' @param detect_sig_thresh Numeric minimum number of comparisons where pathway must be significant (p <= 0.05) to be included (default 0)
#' @param return_list Logical whether to return a summary tibble instead of plot (default FALSE)
#' @param hccol Logical indicating whether to perform hierarchical clustering on columns (default FALSE)
#' @param hcrow Logical indicating whether to perform hierarchical clustering on rows (default FALSE)
#' @param flip_col Logical whether to reverse the column clustering order (default FALSE)
#' @param flip_row Logical whether to reverse the row clustering order (default FALSE)
#' @param show_legend Logical indicating whether to show legends (default TRUE)
#' @param plot_title Optional title for the plot (default NULL)
#' @param save_path Optional file path to save the plot (default NULL)
#' @param plot_width Numeric width for saved plot in inches (default NULL - auto-calculated)
#' @param plot_height Numeric height for saved plot in inches (default NULL - auto-calculated)
#' @param dpi Numeric resolution for saved plot (default 600)
#' @param color_scale Character string for color scheme: "blue" or "red" (default "blue")
#' @param background Background color for saved plot (default "transparent")
#'
#' @return A ggplot2 object or tibble (if return_list = TRUE) with pathway_name, min_p_value, and metabolite columns
#'
#' @examples
#' \dontrun{
#' enrichment_data <- data.frame(
#'   pathway_name = rep(c("Pathway A", "Pathway B"), 2),
#'   Comparisons = rep(c("Group1", "Group2"), each = 2),
#'   enrichment_factor = c(2.5, 3.1, 1.8, 2.9),
#'   p_fisher = c(0.01, 0.02, 0.03, 0.015)
#' )
#'
#' plot <- plot_mummichog_columns(enrichment_data)
#' }
#'
#' @export
plot_mummichog_columns <- function(
    enrichment_data,
    p_display_threshold = 0.05,
    p_color = c(0.1, 0.05, 0.01),
    detect_thresh = 1,
    detect_sig_thresh = 0,
    return_list = FALSE,
    hccol = FALSE,
    hcrow = FALSE,
    flip_col = FALSE,
    flip_row = FALSE,
    show_legend = TRUE,
    plot_title = NULL,
    save_path = NULL,
    plot_width = NULL,
    plot_height = NULL,
    dpi = 600,
    color_scale = "blue",
    background = "transparent") {
  
  # Load required libraries
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(forcats)
  
  # Calculate scale factor: 8/14 ≈ 0.571 to match other plots
  scale_factor <- 8/14
  
  # Source clean_pathway_names_corr function
  if (!exists("clean_pathway_names_corr")) {
    source("R/Utilities/Helpers/clean_pathway_names_corr.R")
  }
  
  # Process data
  enrichment_data <- enrichment_data |>
    # Filter by p-value threshold (p_value is already -log10(p))
    filter(p_fisher >= -log10(p_display_threshold)) |>
    # Clean pathway names with correlation-specific function
    mutate(pathway_name = clean_pathway_names_corr(pathway_name))
  
  # Filter by detection threshold: keep only pathways detected in >= detect_thresh comparisons
  if (detect_thresh > 1) {
    pathway_counts <- enrichment_data |>
      group_by(pathway_name) |>
      summarise(n_detected = n_distinct(Comparisons), .groups = "drop") |>
      filter(n_detected >= detect_thresh)
    
    enrichment_data <- enrichment_data |>
      filter(pathway_name %in% pathway_counts$pathway_name)
  }
  
  # Filter by significant detection threshold: keep only pathways significant in >= detect_sig_thresh comparisons
  if (detect_sig_thresh > 0) {
    sig_threshold <- -log10(0.05)  # Significance at p <= 0.05
    pathway_sig_counts <- enrichment_data |>
      filter(p_fisher >= sig_threshold) |>
      group_by(pathway_name) |>
      summarise(n_sig = n_distinct(Comparisons), .groups = "drop") |>
      filter(n_sig >= detect_sig_thresh)
    
    enrichment_data <- enrichment_data |>
      filter(pathway_name %in% pathway_sig_counts$pathway_name)
  }
  
  # If return_list is TRUE, return summary tibble instead of plot
  if (return_list) {
    pathway_summary <- enrichment_data |>
      group_by(pathway_name) |>
      slice_max(p_fisher, n = 1, with_ties = FALSE) |>
      ungroup() |>
      select(
        pathway_name,
        min_p_value = p_fisher,
        metabolite = Comparisons
      ) |>
      arrange(desc(min_p_value))
    
    return(pathway_summary)
  }
  
  # Hierarchical clustering if requested
  if (hccol || hcrow) {
    # Reshape to wide format for clustering (pathways × comparisons)
    wide_data <- enrichment_data |>
      select(pathway_name, Comparisons, p_fisher) |>
      pivot_wider(names_from = Comparisons, values_from = p_fisher, values_fill = 0)
    
    # Extract matrix for clustering
    mat <- as.matrix(wide_data[, -1])
    rownames(mat) <- wide_data$pathway_name
    
    # Cluster columns (metabolites/comparisons)
    if (hccol && ncol(mat) > 1) {
      col_hc <- hclust(dist(t(mat)), method = "average")
      col_order <- colnames(mat)[col_hc$order]
      # Flip order if requested
      if (flip_col) {
        col_order <- rev(col_order)
      }
      enrichment_data <- enrichment_data |>
        mutate(Comparisons = factor(Comparisons, levels = col_order))
    }
    
    # Cluster rows (pathways)
    if (hcrow && nrow(mat) > 1) {
      row_hc <- hclust(dist(mat), method = "average")
      row_order <- rownames(mat)[row_hc$order]
      # Flip order if requested
      if (flip_row) {
        row_order <- rev(row_order)
      }
      enrichment_data <- enrichment_data |>
        mutate(pathway_name = factor(pathway_name, levels = row_order))
    }
  }
  
  # Ensure factors are set (if not already done by clustering)
  if (!hcrow) {
    enrichment_data <- enrichment_data |>
      mutate(pathway_name = factor(pathway_name, levels = unique(pathway_name)))
  }
  if (!hccol) {
    enrichment_data <- enrichment_data |>
      mutate(Comparisons = factor(Comparisons, levels = unique(Comparisons)))
  }
  
  # Complete all pathway×metabolite combinations with NA for missing values
  # This ensures every facet has data and NA tiles will use na.value color
  enrichment_data <- enrichment_data |>
    complete(pathway_name, Comparisons, fill = list(p_fisher = NA, enrichment_factor = NA))
  
  # Create the plot as heatmap (filled cells instead of dots)
  p <- ggplot(
    enrichment_data,
    aes(x = 0.5, y = 0.5, fill = p_fisher)  # Use fill instead of color
  ) +
    # Filled tiles colored by p-value
    geom_tile(
      width = 1, height = 1,
      color = "grey80", linewidth = 0.3,
      na.rm = TRUE, show.legend = show_legend
    ) +
    facet_grid(
      rows = vars(pathway_name),
      cols = vars(Comparisons),
      switch = "y", drop = FALSE
    ) +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0))
  
  # Fill scale based on color_scale parameter (changed from color to fill)
  if (color_scale == "rb") {
    # Red-to-blue gradient using PGD colors (matches plot_mummichog_enrichment)
    # Maps to -log10(p-value) which is in p_fisher column
    # User specifies p-values in p_color vector: first = darkest blue, last = darkest red
    
    # Transform p_color to -log10 scale
    p_color_log <- -log10(p_color)
    
    # Calculate scale limits - expand slightly to cover all data
    scale_min <- min(c(0, p_color_log)) # Ensure 0 is included for very low significance
    scale_max <- max(p_color_log)
    scale_range <- scale_max - scale_min
    
    # Normalize p_color_log values to 0-1 range for gradient mapping
    gradient_values <- (p_color_log - scale_min) / scale_range
    
    # Generate breaks: 0.5, 1, asterisk, 2+
    asterisk_log <- -log10(0.05)
    all_breaks <- c(0.5, 1, asterisk_log, 2)
    
    # Generate labels: ✱ for p=0.05, "2+" for 2, numbers for others
    color_labels <- sapply(all_breaks, function(x) {
      if (abs(x - asterisk_log) < 0.001) {
        "✱"
      } else if (abs(x - 2) < 0.001) {
        "2+"
      } else {
        as.character(x)
      }
    })
    
    p <- p + scale_fill_gradientn(
      colors = c("#113d6a", "#D8919A", "#800017"), # Blue to red gradient
      values = gradient_values,
      limits = c(scale_min, scale_max),
      breaks = all_breaks,
      labels = color_labels,
      oob = scales::squish,
      na.value = "#113d6a",  # Use dark blue for any NA or out-of-bounds values
      name = "-log10(p-value)\n",
      guide = guide_colorbar(
        reverse = FALSE, # Low (blue) at bottom, high (red) at top
        barheight = unit(10 * scale_factor, "cm"),
        barwidth = unit(0.9 * scale_factor, "cm"),
        ticks.colour = "black",
        ticks.linewidth = 0.5 * scale_factor,
        frame.colour = NA,
        frame.linewidth = 0,
        draw.ulim = TRUE,
        draw.llim = TRUE
      )
    )
  } else if (color_scale == "blue") {
    # Transform p_color to -log10 scale
    p_color_log <- -log10(p_color)
    scale_min <- min(p_color_log)
    scale_max <- max(p_color_log)
    
    # Generate breaks: 0.5, 1, asterisk, 2+
    asterisk_log <- -log10(0.05)
    all_breaks <- c(0.5, 1, asterisk_log, 2)
    color_labels <- sapply(all_breaks, function(x) {
      if (abs(x - asterisk_log) < 0.001) {
        "✱"
      } else if (abs(x - 2) < 0.001) {
        "2+"
      } else {
        as.character(x)
      }
    })
    
    p <- p + scale_fill_gradient(
      low = "#0a2256", high = "#c3dbe9", # Dark (small p) -> light (large p)
      limits = c(scale_min, scale_max),
      breaks = all_breaks,
      labels = color_labels,
      oob = scales::squish,
      na.value = "#0a2256",  # Use dark blue for any NA values
      name = "-log10(p-value)\n",
      guide = guide_colorbar(
        reverse = FALSE,
        barheight = unit(10 * scale_factor, "cm"),
        barwidth = unit(0.9 * scale_factor, "cm"),
        ticks.colour = "black",
        ticks.linewidth = 0.5 * scale_factor,
        frame.colour = NA,
        frame.linewidth = 0,
        draw.ulim = TRUE,
        draw.llim = TRUE
      )
    )
  } else {
    # Red/orange color scheme
    # Transform p_color to -log10 scale
    p_color_log <- -log10(p_color)
    scale_min <- min(p_color_log)
    scale_max <- max(p_color_log)
    
    # Generate breaks: 0.5, 1, asterisk, 2+
    asterisk_log <- -log10(0.05)
    all_breaks <- c(0.5, 1, asterisk_log, 2)
    color_labels <- sapply(all_breaks, function(x) {
      if (abs(x - asterisk_log) < 0.001) {
        "✱"
      } else if (abs(x - 2) < 0.001) {
        "2+"
      } else {
        as.character(x)
      }
    })
    
    p <- p + scale_fill_gradient(
      low = "#801914ff", high = "#e8d5b7ff", # Dark red (small p) -> light (large p)
      limits = c(scale_min, scale_max),
      breaks = all_breaks,
      labels = color_labels,
      oob = scales::squish,
      na.value = "#113d6a",  # Use dark blue for any NA values
      name = "-log10(p-value)\n",
      guide = guide_colorbar(
        reverse = FALSE,
        barheight = unit(10 * scale_factor, "cm"),
        barwidth = unit(0.9 * scale_factor, "cm"),
        ticks.colour = "black",
        ticks.linewidth = 0.5 * scale_factor,
        frame.colour = NA,
        frame.linewidth = 0,
        draw.ulim = TRUE,
        draw.llim = TRUE
      )
    )
  }
  
  # Theme and styling
  p <- p +
    labs(x = NULL, y = NULL, title = plot_title) +
    theme_minimal(base_family = "Arial") +
    theme(
      text = element_text(family = "Arial"),
      panel.grid = element_blank(),
      panel.background = element_blank(),
      plot.background = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.spacing.x = unit(0, "pt"),
      panel.spacing.y = unit(0, "pt"),
      strip.placement = "outside",
      strip.text.x.top = element_text(
        angle = 90, vjust = 0.5, hjust = 0,
        face = "plain", family = "Arial", size = 13 * scale_factor
      ),
      strip.text.y.left = element_text(
        angle = 0, hjust = 1,
        face = "plain", family = "Arial", size = 13 * scale_factor,
        margin = margin(r = 6 * scale_factor)
      ),
      legend.title = element_text(size = 13 * scale_factor, face = "bold", family = "Arial", vjust = 0.3),
      legend.text = element_text(size = 13 * scale_factor, family = "Arial"),
      plot.margin = margin(t = 20 * scale_factor, r = 40 * scale_factor, b = 10 * scale_factor, l = 120 * scale_factor)
    ) +
    coord_cartesian(clip = "off")
  
  # Hide legends if requested
  if (!show_legend) {
    p <- p + theme(legend.position = "none")
  }
  
  # Auto-calculate dimensions if not provided
  if (is.null(plot_width) || is.null(plot_height)) {
    n_comparisons <- length(unique(enrichment_data$Comparisons))
    n_pathways <- length(unique(enrichment_data$pathway_name))
    
    if (is.null(plot_width)) {
      width_per_comparison <- 0.3
      width_base <- if (show_legend) 7 else 4.2
      plot_width <- width_base + (n_comparisons * width_per_comparison)
    }
    
    if (is.null(plot_height)) {
      height_per_pathway <- 0.3
      height_base <- 2
      plot_height <- height_base + (n_pathways * height_per_pathway)
    }
  }
  
  # Save plot if path provided
  if (!is.null(save_path)) {
    # Create directory if it doesn't exist
    save_dir <- dirname(save_path)
    if (!dir.exists(save_dir)) {
      dir.create(save_dir, recursive = TRUE)
    }
    
    ggplot2::ggsave(
      filename = save_path,
      plot = p,
      width = plot_width,
      height = plot_height,
      dpi = dpi,
      units = "in",
      device = ragg::agg_png,
      bg = background
    )
    
    cat("Plot saved to:", save_path, "\n")
    cat("  Dimensions:", plot_width, "x", plot_height, "inches at", dpi, "DPI\n")
  }
  
  return(invisible(p))
}
