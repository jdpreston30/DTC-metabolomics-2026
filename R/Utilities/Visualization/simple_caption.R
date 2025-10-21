#' Simple Text Wrapping and Caption Creation
#'
#' This uses basic text wrapping and multiple text grobs - much simpler!
#'
simple_caption_creation <- function(png_path,
                                   caption_text,
                                   output_path,
                                   caption_size = 11,
                                   figure_width = 8.5,
                                   figure_height = 11,
                                   dpi = 300) {
  
  library(cowplot)
  library(magick)
  library(grid)
  library(stringr)
  
  # Simple text wrapping function
  wrap_text <- function(text, width_chars = 100) {
    # Split into words
    words <- strsplit(text, " ")[[1]]
    lines <- character()
    current_line <- ""
    
    for (word in words) {
      test_line <- if (current_line == "") word else paste(current_line, word)
      
      if (nchar(test_line) > width_chars && current_line != "") {
        lines <- c(lines, current_line)
        current_line <- word
      } else {
        current_line <- test_line
      }
    }
    
    if (current_line != "") {
      lines <- c(lines, current_line)
    }
    
    return(lines)
  }
  
  # Process the caption text to make Figure 1. bold and (A), (B), etc. bold
  process_caption <- function(text) {
    # Make "Figure X." bold
    text <- gsub("^(Figure \\d+\\.)", "**\\1**", text)
    # Make letters in parentheses bold
    text <- gsub("\\(([A-Z])\\)", "**(\\1)**", text)
    return(text)
  }
  
  processed_text <- process_caption(caption_text)
  wrapped_lines <- wrap_text(processed_text, width_chars = 100)
  
  # Read main image
  main_img <- magick::image_read(png_path)
  main_grob <- grid::rasterGrob(as.raster(main_img))
  
  # Calculate proportions
  caption_height_inches <- 1.0
  img_height_prop <- (figure_height - caption_height_inches) / figure_height
  caption_height_prop <- caption_height_inches / figure_height
  
  # Create text grobs for each line
  line_height <- 0.8 / length(wrapped_lines)  # Distribute lines in caption area
  
  # Start building the plot
  combined_plot <- ggdraw() +
    # Add main image
    draw_grob(main_grob, x = 0, y = caption_height_prop, width = 1, height = img_height_prop)
  
  # Add each line of text
  for (i in seq_along(wrapped_lines)) {
    line <- wrapped_lines[i]
    y_pos <- caption_height_prop - (i - 1) * line_height / length(wrapped_lines) - 0.05
    
    # Simple approach: treat ** as bold indicators and create multiple text elements
    if (grepl("\\*\\*", line)) {
      # For now, just remove the ** and make it regular text
      # You could enhance this to actually make parts bold
      clean_line <- gsub("\\*\\*", "", line)
    } else {
      clean_line <- line
    }
    
    combined_plot <- combined_plot +
      draw_label(
        clean_line,
        x = 0.05, y = y_pos,
        hjust = 0, vjust = 1,
        size = caption_size,
        fontfamily = "Arial"
      )
  }
  
  # Save
  ggsave(
    filename = output_path,
    plot = combined_plot,
    width = figure_width,
    height = figure_height,
    dpi = dpi,
    units = "in",
    bg = "white"
  )
  
  cat("Simple caption saved to:", output_path, "\n")
}