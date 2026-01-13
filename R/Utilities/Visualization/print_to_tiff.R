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
