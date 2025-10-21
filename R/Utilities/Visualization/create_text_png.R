#' Create PNG of Formatted Text
#'
#' This function creates a PNG file containing formatted text that can later
#' be read back and used as a grob in cowplot layouts.
#'
#' @param text_content The text content (can include markdown if using ggtext)
#' @param output_path Path for the text PNG file
#' @param text_size Font size
#' @param text_family Font family
#' @param text_color Text color
#' @param bg_color Background color (use "transparent" for transparent)
#' @param width_inches Width of text box in inches
#' @param height_inches Height of text box in inches
#' @param dpi DPI for the PNG
#' @param use_markdown Whether to use markdown formatting (requires ggtext)
#'
create_text_png <- function(text_content,
                           output_path,
                           text_size = 11,
                           text_family = "Arial",
                           text_color = "black",
                           bg_color = "white",
                           width_inches = 8,
                           height_inches = 1.5,
                           dpi = 300,
                           use_markdown = TRUE) {
  
  library(ggplot2)
  library(cowplot)
  
  if (use_markdown) {
    library(ggtext)
    
    # Create plot with markdown support and proper text wrapping
    text_plot <- ggplot() +
      geom_richtext(
        aes(x = 0.05, y = 0.5, label = text_content),
        hjust = 0, vjust = 0.5,  # Left-aligned for better wrapping
        size = text_size / .pt,
        family = text_family,
        color = text_color,
        fill = NA,  # No background fill
        label.color = NA,  # No border
        lineheight = 1.3,
        # Key fix: Use explicit width in inches converted to grid units
        width = unit(width_inches - 0.5, "inches"),  # Leave some margin
        label.padding = unit(c(2, 2, 2, 2), "pt")
      ) +
      xlim(0, 1) + ylim(0, 1) +
      theme_void() +
      theme(
        plot.background = element_rect(fill = bg_color, color = NA),
        panel.background = element_rect(fill = bg_color, color = NA),
        plot.margin = margin(5, 5, 5, 5, "pt")
      )
  } else {
    # For non-markdown, we need to manually wrap text
    library(stringr)
    
    # Simple text wrapping function
    wrap_text_manual <- function(text, width = 100) {
      words <- strsplit(text, " ")[[1]]
      lines <- character()
      current_line <- character()
      
      for (word in words) {
        test_line <- paste(c(current_line, word), collapse = " ")
        if (nchar(test_line) > width && length(current_line) > 0) {
          lines <- c(lines, paste(current_line, collapse = " "))
          current_line <- word
        } else {
          current_line <- c(current_line, word)
        }
      }
      if (length(current_line) > 0) {
        lines <- c(lines, paste(current_line, collapse = " "))
      }
      return(paste(lines, collapse = "\n"))
    }
    
    wrapped_text <- wrap_text_manual(text_content, width = 120)
    
    # Create plot with regular text
    text_plot <- ggplot() +
      geom_text(
        aes(x = 0.05, y = 0.5, label = wrapped_text),
        hjust = 0, vjust = 0.5,
        size = text_size / .pt,
        family = text_family,
        color = text_color,
        lineheight = 1.3
      ) +
      xlim(0, 1) + ylim(0, 1) +
      theme_void() +
      theme(
        plot.background = element_rect(fill = bg_color, color = NA),
        panel.background = element_rect(fill = bg_color, color = NA),
        plot.margin = margin(5, 5, 5, 5, "pt")
      )
  }
  
  # Save the text as PNG
  ggsave(
    filename = output_path,
    plot = text_plot,
    width = width_inches,
    height = height_inches,
    dpi = dpi,
    units = "in",
    bg = bg_color
  )
  
  cat("Text PNG saved to:", output_path, "\n")
  return(output_path)
}

#' Read Text PNG as Grob for Cowplot
#'
#' This function reads a PNG file and converts it to a grob that can be used
#' in cowplot layouts.
#'
#' @param png_path Path to the PNG file
#'
read_text_png_as_grob <- function(png_path) {
  library(magick)
  library(grid)
  
  # Read the PNG
  img <- magick::image_read(png_path)
  
  # Convert to grob
  img_grob <- grid::rasterGrob(as.raster(img))
  
  return(img_grob)
}

# Example usage:
# 
# # Step 1: Create the text PNG
# caption_text <- "**Figure 1.** Your caption with **(A)** bold letters and proper wrapping..."
# create_text_png(
#   text_content = caption_text,
#   output_path = "temp_caption.png",
#   text_size = 11,
#   text_family = "Arial",
#   bg_color = "white",
#   width_inches = 8.5,
#   height_inches = 1.2,
#   use_markdown = TRUE
# )
# 
# # Step 2: Read it back as a grob
# caption_grob <- read_text_png_as_grob("temp_caption.png")
# 
# # Step 3: Use in cowplot
# final_plot <- ggdraw() +
#   draw_grob(your_main_figure_grob, x = 0, y = 0.15, width = 1, height = 0.85) +
#   draw_grob(caption_grob, x = 0, y = 0, width = 1, height = 0.15)