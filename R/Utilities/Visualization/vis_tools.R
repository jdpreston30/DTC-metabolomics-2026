#+ Grid Guide for Positioning (cowplot coordinates)
#' Grid Guide for Positioning (cowplot coordinates)
#' Creates a grid overlay with dashed lines and coordinate labels to assist with
#' precise positioning of plot elements in cowplot figure layouts.
#'
#' @param x_max Maximum x-coordinate for the grid (default: 8.5)
#' @param y_max Maximum y-coordinate for the grid (default: 11)
#' @param interval Spacing between grid lines (default: 0.25)
#' @param label_interval Spacing between coordinate labels (default: 0.5)
#' @param margins Distance from edges for margin lines in inches (default: 0.5)
#'
#' @return List of ggplot2 geom objects for grid overlay
#'
#' @details
#' The function creates a coordinate grid system for precise element positioning:
#' - Red dashed lines at specified intervals
#' - Coordinate labels at specified intervals
#' - Optional black margin lines to show printable area
#' 
#' Coordinates follow cowplot's system where (0,0) is bottom-left.
#'
#' @examples
#' \dontrun{
#'   # Add grid to existing plot
#'   plot + grdgd()
#'   
#'   # Custom grid settings
#'   plot + grdgd(x_max = 10, y_max = 8, interval = 0.5)
#' }
#'
#' @export
grdgd <- function(x_max = 8.5, y_max = 11, interval = 0.25, label_interval = 0.5, margins = 0.5) {
  guide_elements <- list(
    ggplot2::geom_vline(xintercept = seq(0, x_max, interval), color = "red", alpha = 0.3, linetype = "dashed"),
    ggplot2::geom_hline(yintercept = seq(0, y_max, interval), color = "red", alpha = 0.3, linetype = "dashed"),
    # Bottom labels (x-axis)
    ggplot2::annotate("text", x = seq(0, x_max, label_interval), y = 0.2, label = seq(0, x_max, label_interval), size = 3, color = "red"),
    # Left labels (y-axis)
    ggplot2::annotate("text", x = 0.2, y = seq(0, y_max, label_interval), label = seq(0, y_max, label_interval), size = 3, color = "red"),
    # Top labels (x-axis)
    ggplot2::annotate("text", x = seq(0, x_max, label_interval), y = y_max - 0.2, label = seq(0, x_max, label_interval), size = 3, color = "red"),
    # Right labels (y-axis)
    ggplot2::annotate("text", x = x_max - 0.2, y = seq(0, y_max, label_interval), label = seq(0, y_max, label_interval), size = 3, color = "red")
  )
  # Add margin lines if margins is specified
  if (!is.null(margins)) {
    # Add vertical margin lines (left and right)
    guide_elements[[length(guide_elements) + 1]] <- ggplot2::geom_vline(xintercept = margins, color = "black", linewidth = 1)
    guide_elements[[length(guide_elements) + 1]] <- ggplot2::geom_vline(xintercept = x_max - margins, color = "black", linewidth = 1)
    
    # Add horizontal margin lines (bottom and top)
    guide_elements[[length(guide_elements) + 1]] <- ggplot2::geom_hline(yintercept = margins, color = "black", linewidth = 1)
    guide_elements[[length(guide_elements) + 1]] <- ggplot2::geom_hline(yintercept = y_max - margins, color = "black", linewidth = 1)
  }
  return(guide_elements)
}
#+ Figure Labels Generator for Cowplot Layouts
#' Figure Labels Generator for Cowplot Layouts
#' Generates figure panel labels (A, B, C, etc.) at specified coordinates
#' for multi-panel figure layouts using cowplot's coordinate system.
#'
#' @param labels Named list where names are label text and values are coordinate vectors c(x, y)
#' @param size Font size for labels (default: 14)
#' @param fontface Font face for labels (default: "bold")
#' @param fontfamily Font family for labels (default: "Arial")
#' @param hjust Horizontal justification (default: 0)
#'
#' @return List of cowplot::draw_label objects
#'
#' @details
#' Creates publication-ready figure panel labels at precise coordinates.
#' Coordinates use cowplot's system where (0,0) is bottom-left corner.
#' 
#' Common positioning:
#' - Top-left panels: around (0.8, 9.7)
#' - Top-right panels: around (4.3, 9.7)
#' - Bottom panels: adjust y-coordinate accordingly
#'
#' @examples
#' \dontrun{
#'   # Define label positions
#'   labels <- list(
#'     A = c(0.8, 9.7),
#'     B = c(4.3, 9.7),
#'     C = c(0.8, 5.2)
#'   )
#'   
#'   # Add to cowplot layout
#'   final_plot + figure_labels(labels)
#' }
#'
#' @export
figure_labels <- function(labels, size = 14, fontface = "bold", fontfamily = "Arial", hjust = 0) {
  # Convert single label to list format if needed
  if (is.character(labels)) {
    stop("Please provide labels as a named list with x and y coordinates, e.g., list(A = c(0.8, 9.7), B = c(3.7, 9.7))")
  }
  
  # Generate draw_label calls for each label
  label_layers <- list()
  for (name in names(labels)) {
    coords <- labels[[name]]
    label_layers[[length(label_layers) + 1]] <- 
      cowplot::draw_label(name, x = coords[1], y = coords[2], 
                         size = size, fontface = fontface, 
                         fontfamily = fontfamily, hjust = hjust)
  }
  
  return(label_layers)
}
#+ Print plot to PNG with auto-refresh for macOS Preview
#' Print plot to PNG with auto-refresh for macOS Preview
#' @param plot The plot object to print
#' @param filename Name of the PNG file (with or without .png extension)
#' @param width Width in inches (default: 8.5)
#' @param height Height in inches (default: 11)
#' @param dpi Resolution in DPI (default: 600 for high quality)
#' @param output_dir Directory to save the PNG (default: "Figures")
#' @param auto_open Whether to automatically open in Preview on first run (default: TRUE)
#' @param background Background color for the plot (default: "white", can be "transparent")
#' @return Invisible path to the created PNG file
#' @export
print_to_png <- function(plot, filename, width = 8.5, height = 11, dpi = 600,
                          output_dir = "Outputs/Figures", auto_open = TRUE, background = "white") {
  # Ensure filename has .png extension
  if (!grepl("\\.png$", filename, ignore.case = TRUE)) {
    filename <- paste0(filename, ".png")
  }

  # Create full path
  filepath <- file.path(output_dir, filename)

  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Check if file already exists (for auto-open logic)
  file_exists <- file.exists(filepath)

  # Save the plot as PNG
  ggplot2::ggsave(
    filename = filepath,
    plot = plot,
    width = width,
    height = height,
    dpi = dpi,
    units = "in",
    device = "png",
    bg = background
  )

  # Auto-open in Preview only on first run (or if specified)
  if (auto_open && !file_exists) {
    system(paste("open", shQuote(filepath)))
    cat("PNG saved and opened in Preview:", filepath, "\n")
    cat("Preview will auto-refresh when you re-run this function!\n")
  } else {
    cat("PNG updated:", filepath, "\n")
  }

  # Return path invisibly
  invisible(filepath)
}
#+ Print plot to TIFF with auto-refresh for macOS Preview
#' Print plot to TIFF with auto-refresh for macOS Preview
#'
#' @param plot The plot object to print
#' @param filename Name of the TIFF file (with or without .tiff extension)
#' @param width Width in inches (default: 8.5)
#' @param height Height in inches (default: 11)
#' @param dpi Resolution in DPI (default: 600 for high quality)
#' @param output_dir Directory to save the TIFF (default: "Figures")
#' @param auto_open Whether to automatically open in Preview on first run (default: TRUE)
#' @param background Background color for the plot (default: "white", can be "transparent")
#' @return Invisible path to the created TIFF file
#' @export
print_to_tiff <- function(plot, filename, width = 8.5, height = 11, dpi = 600,
                          output_dir = "Outputs/Figures", auto_open = TRUE, background = "white") {
  # Ensure filename has .tiff extension
  if (!grepl("\\.tiff?$", filename, ignore.case = TRUE)) {
    filename <- paste0(filename, ".tiff")
  }

  # Create full path
  filepath <- file.path(output_dir, filename)

  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Check if file already exists (for auto-open logic)
  file_exists <- file.exists(filepath)

  # Save the plot as TIFF with LZW compression
  ggplot2::ggsave(
    filename = filepath,
    plot = plot,
    width = width,
    height = height,
    dpi = dpi,
    units = "in",
    device = ragg::agg_tiff,
    compression = "lzw",
    bg = background
  )

  # Auto-open in Preview only on first run (or if specified)
  if (auto_open && !file_exists) {
    system(paste("open", shQuote(filepath)))
    cat("TIFF saved and opened in Preview:", filepath, "\n")
    cat("Preview will auto-refresh when you re-run this function!\n")
  } else {
    cat("TIFF updated:", filepath, "\n")
  }

  # Return path invisibly
  invisible(filepath)
}
#+ Spacing Helper for Vertical Centering
#' Spacing Helper for Vertical Centering
#' Calculates the new y position needed to center an element when top and bottom
#' spacing is unequal. Useful for aligning legends, labels, and other elements.
#'
#' @param top_space Space from top edge in pixels (or units before DPI conversion)
#' @param bottom_space Space from bottom edge in pixels (or units before DPI conversion)
#' @param current_y Current y position in inches
#' @param dpi DPI conversion factor (default: 300)
#'
#' @return Prints the new centered y position
#'
#' @details
#' When an element has unequal spacing from top and bottom edges, this function
#' calculates the adjustment needed to center it. The total space is divided
#' equally, and the element is repositioned to achieve equal spacing on both sides.
#' 
#' Calculation: 
#' - Total space = top_space + bottom_space
#' - Target spacing = total_space / 2
#' - Adjustment = (target_spacing - current_top_space) / dpi
#' - In cowplot: increasing top_space means moving DOWN (decreasing y)
#'
#' @examples
#' \dontrun{
#'   # Legend has 227 pixels from top, 480 pixels from bottom
#'   # Current y position is 1.888333333
#'   # Total = 707, so each side should be 353.5
#'   sy(227, 480, 1.888333333)
#'   # Output: New y position = 1.46666666
#'   
#'   # With different DPI
#'   sy(150, 300, 5.2, dpi = 150)
#' }
#'
#' @export
sy <- function(top_space, bottom_space, current_y, dpi = 300) {
  # Calculate total space and target equal spacing
  total_space <- top_space + bottom_space
  target_spacing <- total_space / 2
  
  # Calculate adjustment needed (negative = move down, positive = move up)
  adjustment_pixels <- target_spacing - top_space
  adjustment_inches <- -adjustment_pixels / dpi
  
  # Calculate new centered position
  new_y <- current_y + adjustment_inches
  
  # Print result with formatting
  cat("\033[1;31mNew y position = ", format(new_y, nsmall = 9), "\033[0m\n", sep = "")
  
  # Return invisibly
  invisible(new_y)
}