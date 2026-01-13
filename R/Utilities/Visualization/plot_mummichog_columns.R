#' Create Multi-Column Pathway Enrichment Dot Plot from Enrichment Data
#'
#' This function creates a dot plot visualization for multiple pathway enrichment sets
#' arranged in columns. Based on plot_mummichog_enrichment but accepts pre-processed
#' enrichment data with a column indicator for faceting.
#'
#' @param enrichment_data Data frame with columns: pathway_name, Comparisons, enrichment_factor, p_fisher
#' @param p_threshold Numeric p-value threshold for significance (default 0.05)
#' @param enrichment_cap Numeric maximum value to cap enrichment factors (default 5)
#' @param size_range Numeric vector of length 2 for dot size range (default c(5, 10))
#' @param size_breaks Numeric vector for size scale breaks (auto-generated if NULL)
#' @param show_legend Logical indicating whether to show legends (default TRUE)
#' @param plot_title Optional title for the plot (default NULL)
#' @param save_path Optional file path to save the plot (default NULL)
#' @param plot_width Numeric width for saved plot in inches (default NULL - auto-calculated)
#' @param plot_height Numeric height for saved plot in inches (default NULL - auto-calculated)
#' @param dpi Numeric resolution for saved plot (default 600)
#' @param color_scale Character string for color scheme: "blue" or "red" (default "blue")
#' @param background Background color for saved plot (default "transparent")
#'
#' @return A ggplot2 object
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
    p_threshold = 0.05,
    enrichment_cap = 5,
    size_range = c(5, 10),
    size_breaks = NULL,
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
  
  # Process data
  enrichment_data <- enrichment_data |>
    # Apply enrichment factor cap
    mutate(enrichment_factor = pmin(enrichment_factor, enrichment_cap)) |>
    # Ensure pathway_name is factor with proper order
    mutate(pathway_name = factor(pathway_name, levels = unique(pathway_name)))
  
  # Auto-generate size breaks if not provided
  if (is.null(size_breaks)) {
    size_breaks <- c(enrichment_cap, enrichment_cap * 0.6, enrichment_cap * 0.2)
    size_breaks <- size_breaks[size_breaks > 0]
  }
  
  # Create the plot
  p <- ggplot(
    enrichment_data,
    aes(x = 0.5, y = 0.5, size = enrichment_factor, color = p_fisher)
  ) +
    # One dummy tile per facet
    geom_tile(
      data = data.frame(x = 0.5, y = 0.5),
      aes(x = x, y = y),
      width = 1, height = 1,
      fill = "white", colour = "grey80", linewidth = 0.3,
      inherit.aes = FALSE
    ) +
    geom_point(
      alpha = 0.95, shape = 16, stroke = 0,
      na.rm = TRUE, show.legend = show_legend
    ) +
    facet_grid(
      rows = vars(pathway_name),
      cols = vars(Comparisons),
      switch = "y", drop = FALSE
    ) +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    
    # Size scale - keep limits ascending; reverse legend order via guide
    scale_size_continuous(
      range = size_range * scale_factor,
      limits = c(0, enrichment_cap),
      breaks = size_breaks,
      labels = c(paste0(max(size_breaks), "+"), as.character(size_breaks[-1])),
      name = "Enrichment factor",
      guide = guide_legend(reverse = TRUE)
    )
  
  # Color scale based on color_scale parameter
  if (color_scale == "blue") {
    p <- p + scale_color_gradient(
      low = "#0a2256", high = "#c3dbe9", # Dark (small p) -> light (large p)
      limits = c(0.01, p_threshold),
      oob = scales::squish,
      name = "p-value\n",
      guide = guide_colorbar(
        reverse = TRUE, # 0.01 at top, threshold at bottom
        barheight = unit(5 * scale_factor, "cm"),
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
    p <- p + scale_color_gradient(
      low = "#801914ff", high = "#e8d5b7ff", # Dark red (small p) -> light (large p)
      limits = c(0.01, p_threshold),
      oob = scales::squish,
      name = "p-value\n",
      guide = guide_colorbar(
        reverse = TRUE,
        barheight = unit(5 * scale_factor, "cm"),
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
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.spacing.x = unit(0, "pt"),
      panel.spacing.y = unit(0, "pt"),
      strip.placement = "outside",
      strip.text.x.top = element_text(
        angle = 0, vjust = 1,
        face = "bold", family = "Arial", size = 13 * scale_factor
      ),
      strip.text.y.left = element_text(
        angle = 0, hjust = 1,
        face = "bold", family = "Arial", size = 12 * scale_factor,
        margin = margin(r = 6 * scale_factor)
      ),
      legend.title = element_text(size = 11 * scale_factor, face = "bold", family = "Arial", vjust = 0.3),
      legend.text = element_text(size = 11 * scale_factor, family = "Arial"),
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
