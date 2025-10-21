#' Create Figure Caption using R Markdown
#'
#' This function creates a caption by rendering R Markdown to HTML, then to PNG.
#' This is much simpler and handles wrapping and formatting automatically.
#'
create_caption_with_rmarkdown <- function(png_path,
                                         caption_text,
                                         output_path,
                                         caption_size = 11,
                                         figure_width = 8.5,
                                         figure_height = 11,
                                         dpi = 300) {
  
  library(rmarkdown)
  library(webshot2)  # or webshot
  library(cowplot)
  library(magick)
  library(grid)
  
  # Create temporary markdown file
  temp_md <- tempfile(fileext = ".md")
  temp_html <- tempfile(fileext = ".html")
  temp_caption_png <- tempfile(fileext = ".png")
  
  # Write markdown content
  md_content <- paste0(
    "---\n",
    "output: html_document\n",
    "---\n\n",
    '<div style="font-family: Arial; font-size: ', caption_size, 'pt; ',
    'line-height: 1.4; padding: 10px; width: ', (figure_width - 0.5) * 72, 'px;">',
    "\n\n", caption_text, "\n\n</div>"
  )
  
  writeLines(md_content, temp_md)
  
  # Render to HTML
  render(temp_md, output_file = temp_html, quiet = TRUE)
  
  # Convert HTML to PNG
  webshot2::webshot(
    url = temp_html,
    file = temp_caption_png,
    vwidth = figure_width * 72,  # Convert inches to pixels (roughly)
    vheight = 200,  # Adjust as needed
    zoom = dpi / 96  # Scale for DPI
  )
  
  # Read images
  main_img <- magick::image_read(png_path)
  caption_img <- magick::image_read(temp_caption_png)
  
  # Combine images using magick
  combined <- magick::image_append(c(main_img, caption_img), stack = TRUE)
  
  # Save combined image
  magick::image_write(combined, output_path, density = dpi)
  
  # Clean up
  unlink(c(temp_md, temp_html, temp_caption_png))
  
  cat("Figure with caption saved to:", output_path, "\n")
}