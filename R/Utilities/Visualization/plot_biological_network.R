#' Plot Biological Network with Customizable Label Positioning
#'
#' Creates a biological network visualization where pathways are connected by shared metabolites,
#' with options for custom label positioning and node number display for easy customization.
#'
#' @param network_data List containing node_data and edge_data from create_biological_network()
#' @param output_file Path to save the plot (optional)
#' @param title Plot title
#' @param subtitle Plot subtitle
#' @param node_size_range Vector of min and max node sizes
#' @param text_size Size of pathway labels
#' @param max_overlaps Maximum number of overlapping text labels allowed
#' @param plot_width Plot width in inches
#' @param plot_height Plot height in inches
#' @param dpi Plot resolution
#' @param show_legend Whether to show the legend
#' @param seed Random seed for reproducible layouts
#' @param variable_edge_thickness Whether to make edge thickness proportional to shared metabolites
#' @param edge_thickness_range Vector of min and max edge thickness when variable_edge_thickness = TRUE
#' @param layout_algorithm Layout algorithm: "fr", "kk", or "stress"
#' @param layout_iterations Number of iterations for layout algorithm
#' @param pull_outliers Whether to constrain outlier nodes
#' @param max_distance_from_center Maximum distance from center for outlier constraint
#' @param label_position Default label position: "above" or "below"
#' @param show_node_numbers Whether to display numbers on nodes for customization reference
#' @param custom_label_positions Named list for custom positioning (e.g., list(p1 = "below", p3 = "above"))
#' @param labels_below Numeric vector of pathway numbers to position below (e.g., c(4, 6, 9))
#' @param nudge_labels Named list for custom vjust values (e.g., list(p1 = -3, p5 = 3.5))
#' @param nudge_labels_vert Named list for vertical nudging/vjust (e.g., list(p1 = -3, p5 = 2))
#' @param nudge_labels_horiz Named list for horizontal nudging/hjust (e.g., list(p1 = 0.2, p3 = 0.8))
#' @param color_scale Character string for color scheme: "red" (default) or "blue"
#' @param background_color Background color: "transparent" (default) or "white"
#'
#' @examples
#' # Basic usage
#' plot_biological_network(network_data)
#' 
#' # Show node numbers to identify pathways for customization
#' plot_biological_network(network_data, show_node_numbers = TRUE)
#' 
#' # Custom label positioning for specific pathways
#' plot_biological_network(network_data, 
#'                         label_position = "above",
#'                         custom_label_positions = list(p1 = "below", p3 = "below"))
#' 
#' # Bulk specify labels below (easier syntax)
#' plot_biological_network(network_data, labels_below = c(4, 6, 9))
#' 
#' # Fine-tune individual label distances
#' plot_biological_network(network_data, 
#'                         nudge_labels = list(p1 = -3, p5 = 3.5))
#' 
#' # Separate vertical and horizontal control
#' plot_biological_network(network_data,
#'                         nudge_labels_vert = list(p1 = -3, p5 = 2),
#'                         nudge_labels_horiz = list(p1 = 0.2, p3 = 0.8))
#' 
#' # With blue color scheme
#' plot_biological_network(network_data, color_scale = "blue")
#' 
#' # With transparent background (default)
#' plot_biological_network(network_data, background_color = "transparent")
#' 
#' # With white background
#' plot_biological_network(network_data, background_color = "white")
#'
#'  
plot_biological_network <- function(network_data, 
                                    output_file = NULL,
                                    title = "Biological Network",
                                    subtitle = "Pathways connected by shared metabolites",
                                    node_size_range = c(3, 10),
                                    text_size = 3,
                                    max_overlaps = 20,
                                    plot_width = 12,
                                    plot_height = 10,
                                    dpi = 300,
                                    show_legend = TRUE,
                                    seed = 42,
                                    variable_edge_thickness = FALSE,
                                    edge_thickness_range = c(0.3, 2.0),
                                    layout_algorithm = "fr",
                                    layout_iterations = 500,
                                    pull_outliers = TRUE,
                                    max_distance_from_center = NULL,
                                    label_position = "above",
                                    show_node_numbers = FALSE,
                                    custom_label_positions = NULL,
                                    labels_below = NULL,
                                    nudge_labels = NULL,
                                    nudge_labels_vert = NULL,
                                    nudge_labels_horiz = NULL,
                                    color_scale = "red",
                                    background_color = "transparent") {
  
  # Load required libraries
  require(igraph)
  require(ggplot2)
  require(ggraph)
  require(ggrepel)
  
  # Network-specific comprehensive pathway name cleaning function
  source("R/Utilities/Helpers/clean_pathway_names_for_network.R")
  
  # Check if network has edges
  if (nrow(network_data$edge_data) == 0) {
    cat("⚠ No biological edges found - cannot create network plot\n")
    return(NULL)
  }
  
  # Prepare node data for igraph (ensure proper column names)
  node_data_clean <- network_data$node_data
  node_data_clean$name <- node_data_clean$id  # igraph expects 'name' column
  
  # Apply comprehensive network-specific pathway name cleaning
  node_data_clean$pathway_label <- clean_pathway_names_for_network(node_data_clean$label)
  
  # Add node numbers to help with customization
  if (show_node_numbers) {
    cat("\nNode reference numbers for custom label positioning:\n")
    for (i in 1:nrow(node_data_clean)) {
      cat(sprintf("p%d = %s\n", i, node_data_clean$label[i]))
    }
    cat("\nUsage examples:\n")
    cat("- Bulk below: labels_below = c(4, 6, 9)\n")
    cat("- Custom positioning: custom_label_positions = list(p1 = \"below\", p3 = \"above\")\n")
    cat("- Fine-tune distances: nudge_labels = list(p1 = -3, p5 = 3.5)\n")
    cat("- Separate controls: nudge_labels_vert = list(p1 = -3), nudge_labels_horiz = list(p1 = 0.2)\n\n")
  }
  
  # Determine label vjust based on position (default or custom)
  node_data_clean$label_vjust <- ifelse(label_position == "above", -2.75, 2.75)
  node_data_clean$label_hjust <- 0.5  # Default center horizontal alignment
  
  # Apply bulk labels_below if provided (easiest method)
  if (!is.null(labels_below)) {
    for (pathway_num in labels_below) {
      if (pathway_num <= nrow(node_data_clean)) {
        node_data_clean$label_vjust[pathway_num] <- 2.75
      }
    }
  }
  
  # Apply custom label positions if provided (overrides labels_below)
  if (!is.null(custom_label_positions)) {
    for (i in 1:nrow(node_data_clean)) {
      node_num <- paste0("p", i)
      if (node_num %in% names(custom_label_positions)) {
        position <- custom_label_positions[[node_num]]
        node_data_clean$label_vjust[i] <- ifelse(position == "above", -2.75, 2.75)
      }
    }
  }
  
  # Apply specific nudge values if provided (overrides everything else)
  if (!is.null(nudge_labels)) {
    for (i in 1:nrow(node_data_clean)) {
      node_num <- paste0("p", i)
      if (node_num %in% names(nudge_labels)) {
        node_data_clean$label_vjust[i] <- nudge_labels[[node_num]]
      }
    }
  }
  
  # Apply vertical nudging (can be used with or instead of nudge_labels)
  if (!is.null(nudge_labels_vert)) {
    for (i in 1:nrow(node_data_clean)) {
      node_num <- paste0("p", i)
      if (node_num %in% names(nudge_labels_vert)) {
        node_data_clean$label_vjust[i] <- nudge_labels_vert[[node_num]]
      }
    }
  }
  
  # Apply horizontal nudging
  if (!is.null(nudge_labels_horiz)) {
    for (i in 1:nrow(node_data_clean)) {
      node_num <- paste0("p", i)
      if (node_num %in% names(nudge_labels_horiz)) {
        node_data_clean$label_hjust[i] <- nudge_labels_horiz[[node_num]]
      }
    }
  }
  
  # Add node numbers for display if requested
  node_data_clean$node_number <- 1:nrow(node_data_clean)
  
  # Create igraph object from the biological network data
  bio_graph <- graph_from_data_frame(
    d = network_data$edge_data[, c("source", "target", "weight")],
    directed = FALSE,
    vertices = node_data_clean
  )
  
  # Set seed for reproducible layout
  set.seed(seed)
  
  # Create layout with specified algorithm and parameters
  if (layout_algorithm == "fr") {
    if (pull_outliers) {
      # Fruchterman-Reingold with aggressive parameters to pull outliers closer
      bio_network_plot <- ggraph(bio_graph, layout = "fr", 
                                 niter = layout_iterations * 2,  # More iterations
                                 area = vcount(bio_graph)^1.2,   # Much smaller area
                                 repulserad = vcount(bio_graph)^1.4)  # Much less repulsion
    } else {
      # Standard Fruchterman-Reingold
      bio_network_plot <- ggraph(bio_graph, layout = "fr", 
                                 niter = layout_iterations,
                                 area = vcount(bio_graph)^1.5,
                                 repulserad = vcount(bio_graph)^1.8)
    }
  } else if (layout_algorithm == "kk") {
    # Kamada-Kawai layout (naturally more compact)
    bio_network_plot <- ggraph(bio_graph, layout = "kk", 
                               maxiter = layout_iterations)
  } else if (layout_algorithm == "stress") {
    # Stress majorization (often good for outliers)
    bio_network_plot <- ggraph(bio_graph, layout = "stress")
  } else {
    # Default layout
    bio_network_plot <- ggraph(bio_graph, layout = layout_algorithm)
  }
  
  # Apply max distance constraint if specified
  if (!is.null(max_distance_from_center)) {
    # Get the layout coordinates directly from the graph
    layout_coords <- create_layout(bio_graph, layout = layout_algorithm)
    
    # Find center of layout
    center_x <- mean(layout_coords$x, na.rm = TRUE)
    center_y <- mean(layout_coords$y, na.rm = TRUE)
    
    # Calculate distances from center
    distances <- sqrt((layout_coords$x - center_x)^2 + (layout_coords$y - center_y)^2)
    
    # Cap distances at max_distance_from_center
    outliers <- which(distances > max_distance_from_center)
    if (length(outliers) > 0) {
      for (i in outliers) {
        # Calculate direction vector
        direction_x <- layout_coords$x[i] - center_x
        direction_y <- layout_coords$y[i] - center_y
        
        # Normalize and scale to max distance
        current_distance <- sqrt(direction_x^2 + direction_y^2)
        scale_factor <- max_distance_from_center / current_distance
        
        # Update positions
        layout_coords$x[i] <- center_x + direction_x * scale_factor
        layout_coords$y[i] <- center_y + direction_y * scale_factor
      }
      
      # Recreate the plot with constrained layout
      bio_network_plot <- ggraph(layout_coords)
      cat("✓ Capped", length(outliers), "outlier nodes at max distance", max_distance_from_center, "\n")
    } else {
      # No outliers found, use original layout
      bio_network_plot <- ggraph(layout_coords)
    }
  }
  
  # Add the plot layers
  bio_network_plot <- bio_network_plot +
    {if(variable_edge_thickness) {
      geom_edge_link(aes(width = weight), color = "gray60", alpha = 0.6, show.legend = FALSE)
    } else {
      geom_edge_link(color = "gray60", width = 1, alpha = 0.6)
    }} +
    {if(variable_edge_thickness) {
      scale_edge_width_continuous(range = edge_thickness_range, name = "Shared Metabolites")
    }} +
    geom_node_point(aes(size = hits_sig, color = -log10(pvalue)), alpha = 1.0) +
    geom_node_text(aes(label = pathway_label, vjust = label_vjust, hjust = label_hjust), 
                   size = text_size, 
                   fontface = "bold",
                   family = "Arial") +
    {if(show_node_numbers) {
      geom_node_text(aes(label = node_number), 
                     size = text_size + 1, 
                     fontface = "bold",
                     family = "Arial",
                     color = "white",
                     vjust = 0.5,
                     hjust = 0.5)
    }} +
    scale_size_continuous(range = node_size_range, name = "Significant Hits") +
    scale_color_gradientn(
      colors = if(color_scale == "blue") {
        c("#c3dbe9", "#4a90e2", "#0a2256", "#0a2256")
      } else if(color_scale == "rb") {
        c("#113d6a", "#D8919A", "#800017", "#800017")
      } else {
        c("#e8d5b7ff", "#e49c30ff", "#801914ff", "#801914ff")  # Light at 1.0, orange at sig, dark red at 2.0, same dark red to 3.0
      },
      values = c(0, 0.3, 0.75, 1),  # 1.0 to 1.301 to 2.0 (plateau starts) to 2.5 (plateau continues)
      limits = c(1, 2.0),  # Fixed range from 1 to 2.0
      breaks = c(1, 1.301, 1.5, 1.75, 2.0),  # Fixed ticks: 1, asterisk, 1.5, 2.0
      labels = c("1", "✱", "1.5", "1.75", "2.0"),  # Custom labels with centered heavy asterisk at 0.05
      oob = scales::squish,  # Values outside range get clamped to edge colors
      name = "-log10(p-value)\n"
    ) +
    # Expand plot limits to prevent text clipping
    scale_x_continuous(expand = expansion(mult = 0.15)) +
    scale_y_continuous(expand = expansion(mult = 0.15)) +
    theme_graph(base_family = "Arial") +
    labs(title = NULL, subtitle = NULL) +  # Remove title and subtitle entirely
    theme(legend.position = if(show_legend) "right" else "none",
          text = element_text(family = "Arial"),
          plot.title = element_text(family = "Arial", hjust = 0.5, margin = margin(t = 120, b = 15)),  # Much larger top margin for title
          plot.subtitle = element_text(family = "Arial", hjust = 0.5, margin = margin(t = 15, b = 40)),  # Larger margins for subtitle
          legend.title = element_text(family = "Arial"),
          legend.text = element_text(family = "Arial"),
          plot.background = element_rect(fill = background_color, color = NA),  # Set background color
          panel.background = element_rect(fill = background_color, color = NA),  # Set panel background
          plot.margin = margin(160, 20, 20, 20, "pt"))  # Much larger top margin for entire plot
          # plot.margin = margin(20, 20, 20, 20, "pt"))  # Add margins around entire plot
  
  # Save the plot if output file specified
  if (!is.null(output_file)) {
    ggsave(output_file, bio_network_plot, 
           width = plot_width, height = plot_height, dpi = dpi)
    cat("✓ Biological network plot saved to:", output_file, "\n")
  }
  
  # Return the plot object
  return(bio_network_plot)
}
