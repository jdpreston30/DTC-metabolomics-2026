#' Create volcano plot from analysis results
#'
#' @param volcano_results Results object from run_volcano function
#' @param up_color Color for upregulated features (default: "#800017")
#' @param down_color Color for downregulated features (default: "#113d6a")
#' @param x_limits Optional vector of length 2 specifying x-axis limits c(min, max). If NULL, auto-scales.
#' @param y_limits Optional vector of length 2 specifying y-axis limits c(min, max). If NULL, auto-scales.
#' @param point_size Size of the points (default: 0.5)
#' @return ggplot object
#' @export
plot_volcano <- function(volcano_results,
                        up_color = "#800017",
                        down_color = "#113d6a",
                        x_limits = NULL,
                        y_limits = NULL,
                        point_size = 0.5) {
  
  # Extract data from results
  volcano_data <- volcano_results$volcano_data
  fc_threshold <- volcano_results$fc_threshold
  p_threshold <- volcano_results$p_threshold
  up_label <- volcano_results$up_label
  down_label <- volcano_results$down_label



  # Create new consistent labels
  new_up_label <- up_label
  new_down_label <- down_label
  
  # Set legend order: Up, Down, Not Significant
  volcano_data$Legend <- factor(
    volcano_data$Legend,
    levels = c(up_label, down_label, "Not Significant")
  )
  
  # ---- Handle axis limits and warnings ----
  # X-axis limits
  if (!is.null(x_limits)) {
    if (length(x_limits) != 2) {
      stop("x_limits must be a vector of length 2: c(min, max)")
    }
    x_range <- x_limits
    # Check for data outside x-axis range
    x_outside <- volcano_data$log2_fc < x_limits[1] | volcano_data$log2_fc > x_limits[2]
    if (any(x_outside, na.rm = TRUE)) {
      n_outside <- sum(x_outside, na.rm = TRUE)
      warning(paste(n_outside, "data points fall outside specified x-axis range",
                    paste0("(", x_limits[1], ", ", x_limits[2], ")")))
    }
  } else {
    # Auto-scale x-axis
    x_range <- c(min(volcano_data$log2_fc, na.rm = TRUE), 
                 max(volcano_data$log2_fc, na.rm = TRUE))
  }
  
  # Y-axis limits
  if (!is.null(y_limits)) {
    if (length(y_limits) != 2) {
      stop("y_limits must be a vector of length 2: c(min, max)")
    }
    y_range <- y_limits
    # Check for data outside y-axis range
    y_outside <- volcano_data$neg_log10_p < y_limits[1] | volcano_data$neg_log10_p > y_limits[2]
    if (any(y_outside, na.rm = TRUE)) {
      n_outside <- sum(y_outside, na.rm = TRUE)
      warning(paste(n_outside, "data points fall outside specified y-axis range",
                    paste0("(", y_limits[1], ", ", y_limits[2], ")")))
    }
  } else {
    # Auto-scale y-axis
    y_range <- c(-0.25, max(volcano_data$neg_log10_p, na.rm = TRUE) * 1.1)
  }
  
  # ---- Create plot ----
  volcano_plot <- ggplot2::ggplot(
    volcano_data,
    ggplot2::aes(x = log2_fc, y = neg_log10_p, color = Legend)
  ) +
    ggplot2::geom_point(size = point_size, na.rm = TRUE) +
    ggplot2::scale_color_manual(
      values = setNames(c("gray70", up_color, down_color), 
                       c("Not Significant", up_label, down_label)),
      labels = setNames(c("Not Significant", new_up_label, new_down_label),
                       c("Not Significant", up_label, down_label)),
      name = NULL
    ) +
    ggplot2::theme_minimal(base_family = "Arial") +
    ggplot2::labs(
      x = expression(bold("log")[2]*bold("(Fold Change)")),
      y = expression(bold("-log")[10]*bold("(p)"))
    ) +
    ggplot2::theme(
      axis.title.x = ggplot2::element_text(size = 12.5, face = "bold", color = "black"),
      axis.title.y = ggplot2::element_text(size = 12.5, face = "bold", color = "black")
    ) +
    ggplot2::geom_hline(yintercept = -log10(p_threshold), linetype = "dashed", color = "black") +
    ggplot2::geom_vline(xintercept = c(-fc_threshold, fc_threshold), linetype = "dashed", color = "black") +
    ggplot2::scale_x_continuous(
      limits = x_range,
      breaks = function(x) {
        range_val <- diff(range(x))
        if (range_val <= 2) seq(ceiling(x[1]), floor(x[2]), 0.5)
        else if (range_val <= 6) seq(ceiling(x[1]), floor(x[2]), 1)
        else seq(ceiling(x[1]), floor(x[2]), 2)
      },
      minor_breaks = function(x) {
        range_val <- diff(range(x))
        if (range_val <= 2) seq(ceiling(x[1]), floor(x[2]), 0.25)
        else if (range_val <= 6) seq(ceiling(x[1]), floor(x[2]), 0.5)
        else seq(ceiling(x[1]), floor(x[2]), 1)
      },
      expand = c(0, 0)
    ) +
    ggplot2::scale_y_continuous(
      limits = y_range,
      breaks = function(x) {
        range_val <- diff(range(x))
        if (range_val <= 3) seq(0, ceiling(x[2]), 0.5)
        else if (range_val <= 6) seq(0, ceiling(x[2]), 1)
        else seq(0, ceiling(x[2]), 2)
      },
      minor_breaks = function(x) {
        range_val <- diff(range(x))
        if (range_val <= 3) seq(0, ceiling(x[2]), 0.25)
        else if (range_val <= 6) seq(0, ceiling(x[2]), 0.5)
        else seq(0, ceiling(x[2]), 1)
      },
  expand = c(0, 0)
    ) +
    ggplot2::theme(
      # Remove aspect ratio - let user control via export

      
      # Legend styling - top left corner (using robust approach for all ggplot2 versions)
      legend.position = c(0.02, 0.98),  # Position inside plot area
      legend.justification = c(0, 1),   # Anchor at top left of legend box
      legend.direction = "vertical",
      legend.box = "vertical", 
      legend.background = ggplot2::element_blank(),  # Remove box around legend
      legend.box.background = ggplot2::element_blank(),  # Also remove box background
      legend.box.margin = ggplot2::margin(0, 0, 0, 0),
      legend.margin = ggplot2::margin(t = 0, r = 0, b = 0, l = 0),
      legend.key = ggplot2::element_blank(),  # Remove background from legend keys
      legend.key.width = grid::unit(0.3, "cm"),
      legend.key.height = grid::unit(0.3, "cm"), 
      legend.key.size = grid::unit(0.3, "cm"),
      legend.spacing.y = grid::unit(0.05, "cm"),  # Tighter spacing
      legend.text = ggplot2::element_text(size = 7, face = "bold"),  # Match PCA size
      legend.title = ggplot2::element_blank(),
      
      # Make legend positioning more robust against patchwork interference
      legend.box.just = "left",
      
      # Axis styling - matching PCA sizes
      axis.title = ggplot2::element_text(size = 12.5, face = "bold", color = "black"),
      axis.title.x = ggplot2::element_text(size = 12.5, face = "bold", color = "black", margin = ggplot2::margin(t = 5)),
      axis.title.y = ggplot2::element_text(size = 12.5, face = "bold", color = "black", margin = ggplot2::margin(r = 5), hjust = 0.5),
      axis.text = ggplot2::element_text(size = 11, face = "bold", color = "black"),
      
  # Panel styling - clean background, keep only bottom and left axes lines
  panel.grid.major = ggplot2::element_blank(),
  panel.grid.minor = ggplot2::element_blank(), 
  panel.background = ggplot2::element_blank(),
  axis.line.x.bottom = ggplot2::element_line(color = "black", linewidth = 0.6),
  axis.line.y.left = ggplot2::element_line(color = "black", linewidth = 0.6),
  axis.ticks = ggplot2::element_line(color = "black", linewidth = 0.6),
  axis.ticks.length = ggplot2::unit(0.15, "cm"),
    ) +
    ggplot2::guides(color = ggplot2::guide_legend(
      override.aes = list(shape = 16, size = 0.75) # legend dots 1.5x larger than plot dots
    ))

  return(volcano_plot)
}
