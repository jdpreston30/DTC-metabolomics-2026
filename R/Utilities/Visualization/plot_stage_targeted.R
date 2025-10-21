#' Create Individual Feature Bar Plots for Stage Analysis
#'
#' This function creates publication-style bar plots for individual metabolomic features,
#' showing stage comparisons with means, individual data points, and statistical results.
#' Specifically designed for stage analysis (Stage I/II vs Stage III/IV).
#'
#' @param feature_table Data frame with stage_bin column and feature columns (samples as rows, features as columns)
#' @param metadata_table Tibble with columns: display_name (or Identified Name), p_value, p_value_fdr, feature
#' @param base_family Font family for plots (default: "Arial")
#' @param include_individual_points Whether to show individual data points (default: TRUE)
#' @param jitter_width Width of jitter for individual points (default: 0.15)
#' @param undo_log Whether to undo log2 transformation (2^x) for display (default: FALSE)
#' @param text_scale Scaling factor for all text elements (default: 1.0, use 0.8 for 80% size, etc.)
#' @param use_identified_name Whether to use Identified Name as title instead of feature name (default: TRUE)
#'
#' @return List of ggplot objects, one for each feature
#'
#' @examples
#' \dontrun{
#'   # Create individual plots
#'   plots <- plot_stage_targeted(
#'     feature_table = TFT_annot,
#'     metadata_table = annot_results_clean,
#'     include_individual_points = TRUE,
#'     undo_log = TRUE,
#'     text_scale = 0.6
#'   )
#'   
#'   # View individual plot
#'   print(plots[[1]])
#' }
#'
#' @export
plot_stage_targeted <- function(feature_table,
                               metadata_table,
                               base_family = "Arial",
                               include_individual_points = TRUE,
                               jitter_width = 0.15,
                               undo_log = FALSE,
                               text_scale = 1.0,
                               use_identified_name = TRUE) {
  
  # Load required libraries
  library(dplyr)
  library(ggplot2)
  library(purrr)
  library(stringr)
  library(labeling)
  
  # Set dynamic y-label based on undo_log parameter
  y_label <- if (undo_log) {
    "Spectral Intensity"
  } else {
    expression(Log[2]~Intensity)
  }
  
  # Color mapping: Stage I/II = blue, Stage III/IV = red (using your specified colors)
  group_colors <- c("Stage I/II" = "#113d6a", "Stage III/IV" = "#800017")
  group_colors_light <- c("Stage I/II" = "#87A6C7", "Stage III/IV" = "#D8919A")
  
  # Display labels: shorter versions for x-axis
  display_labels <- c("Stage I/II" = "I/II", "Stage III/IV" = "III/IV")
  x_axis_title <- NULL  # Remove x-axis title
  
  # Publication-style theme
  theme_pub_barplot <- function() {
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
  
  # Function to create a single plot
  create_single_plot <- function(feature_name) {
    # Get metadata for this feature
    feature_meta <- metadata_table %>%
      filter(feature == feature_name)
    
    if (nrow(feature_meta) == 0) {
      warning("No metadata found for feature: ", feature_name)
      return(NULL)
    }
    
    # Get feature values and group info
    feature_values <- feature_table[[feature_name]]
    group_vector <- feature_table$stage_bin
    
    # Convert group values to match our color scheme
    group_vector <- ifelse(group_vector == "Early", "Stage I/II", 
                          ifelse(group_vector == "Late", "Stage III/IV", group_vector))
    
    # Remove missing values
    complete_cases <- !is.na(feature_values) & !is.na(group_vector)
    feature_clean <- feature_values[complete_cases]
    group_clean <- group_vector[complete_cases]
    
    # Convert group to factor with correct order (I/II first, III/IV second)
    if (!is.factor(group_clean)) {
      group_clean <- factor(group_clean, levels = c("Stage I/II", "Stage III/IV"))
    }
    
    # Transform data if requested
    if (undo_log) {
      feature_clean <- 2^feature_clean
    }
    
    # Create plot data
    plot_data <- data.frame(
      Group = group_clean,
      Feature_Value = feature_clean
    )
    
    # Calculate summary statistics
    summary_data <- plot_data %>%
      group_by(Group) %>%
      summarise(
        mean_value = mean(Feature_Value, na.rm = TRUE),
        max_value = max(Feature_Value, na.rm = TRUE),
        n_samples = n(),
        .groups = "drop"
      )
    
    # Get p-values and display name
    p_value <- feature_meta$p_value[1]
    p_value_fdr <- feature_meta$p_value_fdr[1]
    
    # Determine plot title
    if (use_identified_name && "Identified Name" %in% colnames(feature_meta) && !is.na(feature_meta$`Identified Name`[1])) {
      plot_title <- feature_meta$`Identified Name`[1]
    } else if ("display_name" %in% colnames(feature_meta) && !is.na(feature_meta$display_name[1])) {
      plot_title <- feature_meta$display_name[1]
    } else {
      plot_title <- feature_name  # fallback to feature name
    }
    
    # Convert daggers to superscript format internally
    if (!is.null(plot_title) && is.character(plot_title)) {
      plot_title <- gsub("†", '*"†"', plot_title)
    }
    
    # Determine y limits with consistent 4-tick structure
    y_max <- max(summary_data$max_value, na.rm = TRUE)
    
    # Calculate a nice round number for the top tick that's above y_max
    # Use powers of 10 approach for clean scientific notation
    top_tick_magnitude <- ceiling(log10(y_max))
    top_tick_base <- ceiling(y_max / (10^(top_tick_magnitude - 1))) * (10^(top_tick_magnitude - 1))
    
    # If that's not enough, try the next increment
    if (top_tick_base < y_max * 1.1) {
      top_tick_base <- top_tick_base + (10^(top_tick_magnitude - 1))
    }
    
    # Set up 4 evenly spaced ticks: 0, 1/3, 2/3, top
    tick_spacing <- top_tick_base / 3
    y_breaks <- c(0, tick_spacing, 2 * tick_spacing, top_tick_base)
    
    # Set y_limit to 1.02 times the top tick value
    y_limit <- top_tick_base * 1.017
    
    # Create the base plot
    p <- ggplot() +
      # Add bars
      geom_col(
        data = summary_data,
        aes(x = Group, y = mean_value, fill = Group, color = Group),
        width = 0.7, 
        linewidth = 0.6,
        alpha = 0.8, 
        na.rm = TRUE
      )
    
    # Add individual data points if requested
    if (include_individual_points) {
      # Create color vectors for points (dark colors) - Stage I/II = blue, Stage III/IV = red
      point_colors <- ifelse(plot_data$Group == "Stage I/II", "#113d6a", "#800017")
      
      p <- p + geom_jitter(
        data = plot_data,
        aes(x = Group, y = Feature_Value),
        color = point_colors,
        width = jitter_width, 
        size = 0.6, 
        alpha = 1, 
        shape = 16,  # Simple solid circles
        show.legend = FALSE
      )
    }
    
    # Apply colors and styling
    p <- p +
      scale_fill_manual(values = group_colors_light, drop = FALSE) +  # Light colors for bars
      scale_color_manual(values = group_colors, drop = FALSE) +
      scale_x_discrete(labels = display_labels) +
      scale_y_continuous(
        expand = expansion(mult = c(0, 0)), 
        limits = c(0, y_limit),
        breaks = y_breaks,
        labels = function(x) {
          # Custom scientific notation: "5e6" instead of "5e+06"
          sapply(x, function(val) {
            if (val == 0) return("0")
            sci_notation <- sprintf("%.0e", val)
            # Remove the "+" and leading zeros from exponent
            cleaned <- gsub("e\\+0*", "E", sci_notation)
            cleaned <- gsub("E0+$", "", cleaned)  # Remove "e0" for exponent 0
            return(cleaned)
          })
        }
      ) +
      theme_pub_barplot() +
      labs(
        x = x_axis_title,
        y = y_label,
        title = {
          # Check if title contains superscript markup
          if (!is.null(plot_title) && is.character(plot_title) && grepl("\\*", plot_title, fixed = TRUE)) {
            tryCatch({
              parse(text = plot_title)
            }, error = function(e) {
              warning("Failed to parse title: ", plot_title, ". Using regular title.")
              plot_title
            })
          } else {
            plot_title
          }
        }
      )
    
    # Add p-value annotations in top right corner
    p_annotation_text <- paste0(
      "p = ", ifelse(p_value < 0.001, "< 0.001", formatC(p_value, format = "f", digits = 3)), ", ",
      "q = ", ifelse(p_value_fdr < 0.001, "< 0.001", formatC(p_value_fdr, format = "f", digits = 3))
    )
    
    # Position in top right corner without expanding plot area
    x_pos <- 2.5  # Slightly right of the second bar but within reasonable bounds
    y_pos <- y_limit * 0.9999
    
    p <- p + annotate("text", 
                     x = x_pos, y = y_pos, 
                     label = p_annotation_text,
                     hjust = 1, vjust = 1,  # Right-align and top-align
                     size = 3 * text_scale,  # Scale the p-value annotation text too
                     fontface = "italic",
                     color = "black")
    
    return(p)
  }
  
  # Get feature columns (all except stage_bin)
  feature_columns <- setdiff(colnames(feature_table), "stage_bin")
  
  # Filter to only features that exist in metadata
  available_features <- intersect(feature_columns, metadata_table$feature)
  
  if (length(available_features) == 0) {
    stop("No matching features found between feature_table and metadata_table")
  }
  
  cat("Creating", length(available_features), "feature plots...\n")
  
  # Create plots for all available features
  plots <- map(available_features, create_single_plot)
  
  # Name the plots using Identified Name if available
  plot_names <- map_chr(available_features, function(feat) {
    meta_row <- metadata_table %>% filter(feature == feat)
    if (nrow(meta_row) > 0) {
      # Use Identified Name if available, otherwise fall back to feature name
      if (use_identified_name && "Identified Name" %in% colnames(metadata_table) && !is.na(meta_row$`Identified Name`[1])) {
        return(meta_row$`Identified Name`[1])
      } else if ("display_name" %in% colnames(metadata_table) && !is.na(meta_row$display_name[1])) {
        return(meta_row$display_name[1])
      } else {
        return(feat)
      }
    } else {
      return(feat)
    }
  })
  
  names(plots) <- plot_names
  
  # Remove any NULL plots
  plots <- plots[!map_lgl(plots, is.null)]
  
  return(plots)
}
