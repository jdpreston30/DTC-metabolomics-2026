smart_preview <- function(plot, width = 8.5, height = 11, scale_factor = 0.6) {
  # Scale down the display size for laptop screens
  display_width = width * scale_factor
  display_height = height * scale_factor
  
  if (length(dev.list()) == 0 || !"quartz" %in% names(dev.list())) {
    quartz(width = display_width, height = display_height, title = "Live Preview")
  }
  print(plot)
}
