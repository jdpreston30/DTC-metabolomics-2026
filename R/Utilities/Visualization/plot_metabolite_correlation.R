#' Plot Metabolite-Metabolite Correlation Scatter Plot
#'
#' Creates a scatter plot showing the correlation between two metabolites,
#' colored by stage, with correlation statistics annotated.
#'
#' @param y_metabolite Display name of the metabolite for y-axis
#' @param x_metabolite Display name of the metabolite for x-axis
#' @param feature_table Data frame with stage_bin column and feature columns (samples as rows)
#' @param metadata_table Tibble with columns: feature, display_name for lookup
#' @param base_family Font family for plots (default: "Arial")
#' @param point_size Size of the points (default: 1.2)
#' @param line_width Width of the correlation line (default: 0.6)
#' @param text_scale Scaling factor for all text elements (default: 1.0)
#' @param show_ci Whether to show confidence interval ribbon (default: FALSE)
#' @param untransform Whether to reverse log2 transformation (2^x) before plotting (default: FALSE)
#' @param minx Minimum value (first tick) on x-axis (default: NULL for auto)
#' @param maxx Maximum value (last tick) on x-axis (default: NULL for auto)
#' @param miny Minimum value (first tick) on y-axis (default: NULL for auto)
#' @param maxy Maximum value (last tick) on y-axis (default: NULL for auto)
#' @param tickx Interval between x-axis ticks (default: NULL for auto)
#' @param ticky Interval between y-axis ticks (default: NULL for auto)
#'
#' @return ggplot object
#'
#' @examples
#' \dontrun{
#'   p <- plot_metabolite_correlation(
#'     y_metabolite = "GMP",
#'     x_metabolite = "Ribose 5-phosphate",
#'     feature_table = TFT_annot_transformed,
#'     metadata_table = direct_corr_metabs
#'   )
#'   print(p)
#' }
#'
#' @export
plot_metabolite_correlation <- function(y_metabolite,
                                       x_metabolite,
                                       feature_table,
                                       metadata_table,
                                       base_family = "Arial",
                                       point_size = 1.2,
                                       line_width = 0.6,
                                       text_scale = 1.0,
                                       untransform = FALSE,
                                       show_ci = FALSE,
                                       minx = NULL,
                                       maxx = NULL,
                                       miny = NULL,
                                       maxy = NULL,
                                       tickx = NULL,
                                       ticky = NULL) {
  
  # Load required libraries  # Load required libraries
  library(dplyr)
  library(ggplot2)
  
  # Look up feature names from display names
  y_feature <- metadata_table |>
    filter(display_name == y_metabolite) |>
    pull(feature)
  
  x_feature <- metadata_table |>
    filter(display_name == x_metabolite) |>
    pull(feature)
  
  # Check if features were found
  if (length(y_feature) == 0) {
    stop(paste("Y metabolite not found:", y_metabolite))
  }
  if (length(x_feature) == 0) {
    stop(paste("X metabolite not found:", x_metabolite))
  }
  
  # Handle multiple matches (take first)
  if (length(y_feature) > 1) {
    warning(paste("Multiple matches for", y_metabolite, "- using first match"))
    y_feature <- y_feature[1]
  }
  if (length(x_feature) > 1) {
    warning(paste("Multiple matches for", x_metabolite, "- using first match"))
    x_feature <- x_feature[1]
  }
  
  # Extract data
  plot_data <- feature_table |>
    select(stage_bin, x = all_of(x_feature), y = all_of(y_feature)) |>
    filter(!is.na(x), !is.na(y), !is.na(stage_bin))
  
  # Untransform data if requested (reverse log2 transformation)
  if (untransform) {
    plot_data$x <- 2^plot_data$x
    plot_data$y <- 2^plot_data$y
  }
  
  # Print data ranges for user reference
  cat(sprintf("\n--- Axis Ranges for %s vs %s ---\n", y_metabolite, x_metabolite))
  cat(sprintf("X-axis (%s): min = %.2f, max = %.2f\n", x_metabolite, min(plot_data$x, na.rm = TRUE), max(plot_data$x, na.rm = TRUE)))
  cat(sprintf("Y-axis (%s): min = %.2f, max = %.2f\n", y_metabolite, min(plot_data$y, na.rm = TRUE), max(plot_data$y, na.rm = TRUE)))
  cat("---------------------------------------\n\n")
  
  # Validate axis parameters if provided (after printing ranges so user can see actual data)
  if (!is.null(minx) && !is.null(maxx) && !is.null(tickx)) {
    if (maxx <= minx) {
      stop("maxx must be greater than minx")
    }
    if (tickx <= 0) {
      stop("tickx must be positive")
    }
    x_range <- maxx - minx
    n_ticks <- x_range / tickx
    if (abs(n_ticks - round(n_ticks)) > 1e-10) {
      stop(sprintf("X-axis tick interval (%s) does not evenly divide the range from minx (%s) to maxx (%s). Range = %s.", 
                   tickx, minx, maxx, x_range))
    }
  }
  
  if (!is.null(miny) && !is.null(maxy) && !is.null(ticky)) {
    if (maxy <= miny) {
      stop("maxy must be greater than miny")
    }
    if (ticky <= 0) {
      stop("ticky must be positive")
    }
    y_range <- maxy - miny
    n_ticks <- y_range / ticky
    if (abs(n_ticks - round(n_ticks)) > 1e-10) {
      stop(sprintf("Y-axis tick interval (%s) does not evenly divide the range from miny (%s) to maxy (%s). Range = %s.", 
                   ticky, miny, maxy, y_range))
    }
  }
  
  # Color mapping: Stage I/II = dark blue, Stage III/IV = dark red
  stage_colors <- c("Stage I/II" = "#113d6a", "Stage III/IV" = "#800017")
  
  # Calculate correlation on WHOLE dataset
  cor_test <- cor.test(plot_data$x, plot_data$y, method = "pearson")
  r_value <- cor_test$estimate
  p_value <- cor_test$p.value
  n <- nrow(plot_data)
  
  # Fit linear model for whole dataset
  lm_fit <- lm(y ~ x, data = plot_data)
  
  # Create prediction data for line
  x_range <- range(plot_data$x, na.rm = TRUE)
  pred_data <- data.frame(x = seq(x_range[1], x_range[2], length.out = 100))
  pred_data$y <- predict(lm_fit, newdata = pred_data)
  
  # For confidence interval if requested
  if (show_ci) {
    pred_ci <- predict(lm_fit, newdata = pred_data, interval = "confidence", level = 0.95)
    pred_data$ymin <- pred_ci[, "lwr"]
    pred_data$ymax <- pred_ci[, "upr"]
  }
  
  # Publication-style theme (matching plot_stage_targeted)
  theme_pub_scatter <- function() {
    theme_minimal(base_family = base_family) +
      theme(
        # Panel styling - clean background, keep only bottom and left axes lines
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        plot.background = element_blank(),
        axis.line.x.bottom = element_line(color = "black", linewidth = 0.6),
        axis.line.y.left = element_line(color = "black", linewidth = 0.6),
        
        # Axis styling (scaled text)
        axis.ticks = element_line(color = "black", linewidth = 0.6),
        axis.ticks.length = unit(0.15, "cm"),
        axis.text = element_text(size = 11 * text_scale, face = "bold", color = "black"),
        axis.text.x = element_text(size = 11 * text_scale, face = "bold", color = "black", angle = 0, vjust = 0.5, hjust = 0.5),
        axis.text.y = element_text(size = 11 * text_scale, face = "bold", color = "black"),
        axis.title = element_text(size = 12 * text_scale, face = "bold", color = "black"),
        axis.title.x = element_text(size = 12 * text_scale, face = "bold", color = "black", margin = margin(t = 5)),
        axis.title.y = element_text(size = 12 * text_scale, face = "bold", color = "black", margin = margin(r = 5), hjust = 0.5),
        
        # Plot title (scaled)
        plot.title = element_text(size = 12 * text_scale, face = "bold", hjust = 0.5, color = "black"),
        
        # Legend
        legend.position = "none"
      )
  }
  
  # Create base plot
  p <- ggplot(plot_data, aes(x = x, y = y)) +
    # Add confidence interval if requested
    {if (show_ci) geom_ribbon(
      data = pred_data,
      aes(x = x, y = y, ymin = ymin, ymax = ymax),
      fill = "gray70",
      alpha = 0.3,
      inherit.aes = FALSE
    )} +
    # Add regression line (for whole dataset, black)
    geom_line(
      data = pred_data,
      aes(x = x, y = y),
      color = "black",
      linewidth = line_width,
      inherit.aes = FALSE
    ) +
    # Add points colored by stage
    geom_point(
      aes(color = stage_bin),
      size = point_size,
      alpha = 0.8,
      shape = 16
    ) +
    # Color scale
    scale_color_manual(
      values = stage_colors,
      name = NULL,
      labels = c("Stage I/II" = "I/II", "Stage III/IV" = "III/IV")
    ) +
    # Force integer scales on both axes
    scale_x_continuous(
      labels = function(x) as.character(as.integer(round(x))),
      expand = expansion(mult = c(0, 0.02)),
      limits = if (!is.null(minx) && !is.null(maxx)) c(minx, maxx) else NULL,
      breaks = if (!is.null(minx) && !is.null(maxx) && !is.null(tickx)) {
        seq(minx, maxx, by = tickx)
      } else {
        waiver()
      }
    ) +
    scale_y_continuous(
      labels = function(y) as.character(as.integer(round(y))),
      expand = expansion(mult = c(0, 0.02)),
      limits = if (!is.null(miny) && !is.null(maxy)) c(miny, maxy) else NULL,
      breaks = if (!is.null(miny) && !is.null(maxy) && !is.null(ticky)) {
        seq(miny, maxy, by = ticky)
      } else {
        waiver()
      }
    ) +
    # Axis labels
    labs(
      x = x_metabolite,
      y = y_metabolite
    ) +
    # Apply theme
    theme_pub_scatter()
  
  return(p)
}
