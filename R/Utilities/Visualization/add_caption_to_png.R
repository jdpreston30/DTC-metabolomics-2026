#' Add Figure Caption using Simple Text Wrapping
#'
#' Uses R's built-in strwrap() function for reliable text wrapping
#'
add_caption_to_png <- function(png_path, 
                              caption_text, 
                              output_path,
                              caption_size = 11,
                              caption_family = "Arial",
                              figure_width = 8.5,
                              figure_height = 11,
                              dpi = 300,
                              caption_height_inches = 1.0) {
  
  library(cowplot)
  library(ggplot2)
  library(magick)
  library(grid)
  
  # Use R's built-in text wrapping
  wrapped_lines <- strwrap(caption_text, width = 100, simplify = FALSE)[[1]]
  
  # Read main image
  main_img <- magick::image_read(png_path)
  main_grob <- grid::rasterGrob(as.raster(main_img))
  
  # Calculate proportions
  img_height_prop <- (figure_height - caption_height_inches) / figure_height
  caption_height_prop <- caption_height_inches / figure_height
  
  # Create the combined plot
  combined_plot <- ggdraw() +
    # Add main image
    draw_grob(main_grob, x = 0, y = caption_height_prop, width = 1, height = img_height_prop)
  
  # Add each line of wrapped text
  line_spacing <- 0.15  # Increase spacing between lines
  start_y <- caption_height_prop - 0.05  # Start closer to the top of caption area
  
  for (i in seq_along(wrapped_lines)) {
    line_text <- wrapped_lines[i]
    y_position <- start_y - (i - 1) * line_spacing
    
    cat("Adding line", i, "at y =", y_position, ":", substr(line_text, 1, 50), "...\n")
    
    # Make Figure X. bold by detecting it
    if (i == 1 && grepl("^Figure \\d+\\.", line_text)) {
      # Split the first line into bold and regular parts
      figure_part <- regmatches(line_text, regexpr("^Figure \\d+\\.", line_text))
      rest_part <- sub("^Figure \\d+\\.\\s*", "", line_text)
      
      # Add bold "Figure X." part
      combined_plot <- combined_plot +
        draw_label(
          figure_part,
          x = 0.05, y = y_position,
          hjust = 0, vjust = 1,
          size = caption_size,
          fontfamily = caption_family,
          fontface = "bold",
          color = "red"  # Temporarily make it red so we can see it
        )
      
      # Calculate approximate width of bold part to position regular text
      bold_width <- nchar(figure_part) * 0.01  # Increase this estimate
      
      # Add regular text part
      if (rest_part != "") {
        combined_plot <- combined_plot +
          draw_label(
            rest_part,
            x = 0.05 + bold_width, y = y_position,
            hjust = 0, vjust = 1,
            size = caption_size,
            fontfamily = caption_family,
            fontface = "plain",
            color = "blue"  # Temporarily make it blue so we can see it
          )
      }
    } else {
      # Regular line
      combined_plot <- combined_plot +
        draw_label(
          line_text,
          x = 0.05, y = y_position,
          hjust = 0, vjust = 1,
          size = caption_size,
          fontfamily = caption_family,
          color = "black"
        )
    }
    }
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
  
  cat("Figure with wrapped caption saved to:", output_path, "\n")
}

# Example usage:
# add_caption_to_png(
#   png_path = "Figures/fig1.png",
#   caption_text = "Figure 1. Your caption text here...",
#   output_path = "Figures/fig1_with_caption.png",
#   figure_width = 8.5,
#   figure_height = 11,
#   dpi = 300
# )