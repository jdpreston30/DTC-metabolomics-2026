list_tree <- function(path = ".", prefix = "") {
  items <- list.files(path, full.names = TRUE)
  for (i in seq_along(items)) {
    item <- items[i]
    name <- basename(item)
    cat(prefix, if (i == length(items)) "└── " else "├── ", name, "\n", sep = "")
    if (dir.exists(item)) {
      list_tree(item, paste0(prefix, if (i == length(items)) "    " else "│   "))
    }
  }
}