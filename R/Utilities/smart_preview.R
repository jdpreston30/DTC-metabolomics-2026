smart_preview <- function(plot, width = 8.5, height = 11) {
  if (length(dev.list()) == 0 || !"quartz" %in% names(dev.list())) {
    quartz(width = width, height = height, title = "Live Preview")
  }
  print(plot)
}