run_scripts <- function(spec = "06", raw_path) {
  purrr::walk(
    list.files(
      here::here("R", "Utilities"),
      pattern = "\\.[rR]$",
      full.names = TRUE,
      recursive = TRUE
    ),
    source
  )

  script_dir <- here::here("R", "Scripts")
  scripts <- gtools::mixedsort(list.files(script_dir, pattern = "\\.[rR]$", full.names = TRUE))

  # Put checkpoints inside raw_path
  chk_dir <- file.path(raw_path, "Checkpoints")
  if (!dir.exists(chk_dir)) dir.create(chk_dir, recursive = TRUE)

  # Helper: extract prefix (e.g. "01", "02")
  script_ids <- substr(basename(scripts), 1, 2)

  # Parse spec argument
  to_run <- c()
  parts <- unlist(strsplit(spec, ","))
  for (p in parts) {
    if (grepl("-", p)) {
      range <- unlist(strsplit(p, "-"))
      start <- which(script_ids == range[1])
      end <- which(script_ids == range[2])
      if (length(start) == 0 | length(end) == 0) stop("Invalid range: ", p)
      to_run <- c(to_run, seq(start, end))
    } else {
      idx <- which(script_ids == p)
      if (length(idx) == 0) stop("Script not found for ID: ", p)
      to_run <- c(to_run, idx)
    }
  }
  to_run <- unique(sort(to_run))

  # Run scripts
  for (i in to_run) {
    script <- scripts[i]
    script_name <- tools::file_path_sans_ext(basename(script))
    message(">>> Running: ", basename(script))

    tryCatch(
      {
        source(script, local = .GlobalEnv, echo = TRUE, keep.source = TRUE)

        # Save checkpoint after each script
        chk_file <- file.path(chk_dir, paste0(script_name, ".rds"))
        saveRDS(as.list(.GlobalEnv), chk_file)
        message(">>> Checkpoint saved: ", chk_file)
      },
      error = function(e) {
        message("❌ Error in ", basename(script), ": ", e$message)
        stop(e)
      }
    )
  }

  message("✅ Finished running scripts: ", spec)
}
