#' Load all checkpoint .rds files under raw_path/Checkpoints
#' @param raw_path Path containing the "Checkpoints" directory
#' @param envir Environment to load into (default = .GlobalEnv)
load_checkpoints <- function(raw_path, envir = .GlobalEnv) {
  checkpoint_dir <- file.path(raw_path, "Checkpoints")

  if (!dir.exists(checkpoint_dir)) {
    stop("❌ Checkpoint directory not found: ", checkpoint_dir)
  }

  # find all .rds files recursively
  rds_files <- list.files(
    checkpoint_dir,
    pattern = "\\.rds$",
    full.names = TRUE,
    recursive = TRUE
  )

  if (length(rds_files) == 0) {
    stop("❌ No .rds files found in ", checkpoint_dir)
  }

  loaded_names <- c()

  for (f in rds_files) {
    message(">>> Loading checkpoint: ", f)
    objs <- readRDS(f)
    list2env(objs, envir = envir)
    loaded_names <- c(loaded_names, names(objs))
  }

  loaded_names <- unique(loaded_names)
  message(
    "✅ Loaded ", length(loaded_names), " unique objects from ",
    length(rds_files), " checkpoint file(s)"
  )

  invisible(loaded_names)
}
