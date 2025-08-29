load_with_checkpoints <- function(raw_path, files_csv, files_xlsx, checkpoint_dir = here::here("R", "Checkpoints")) {
  # Ensure checkpoint directory exists
  if (!dir.exists(checkpoint_dir)) dir.create(checkpoint_dir, recursive = TRUE)
  rds_file <- file.path(checkpoint_dir, "all_objects.rds")

  if (file.exists(rds_file)) {
    message(">>> Loading from checkpoint: ", rds_file)
    objs <- readRDS(rds_file)
    list2env(objs, envir = .GlobalEnv)
    rm(objs)
  } else {
    message(">>> No checkpoint found, reading raw data...")

    # --- Load CSVs ---
    csv_objs <- purrr::map2(
      files_csv$name, files_csv$file,
      ~ readr::read_csv(file.path(raw_path, paste0(.y, ".csv")))
    )
    names(csv_objs) <- files_csv$name

    # --- Load XLSX ---
    xlsx_objs <- purrr::pmap(
      list(files_xlsx$name, files_xlsx$file, files_xlsx$sheet),
      ~ readxl::read_excel(file.path(raw_path, paste0(..2, ".xlsx")), sheet = ..3)
    )
    names(xlsx_objs) <- files_xlsx$name

    # Combine all objects
    all_objs <- c(csv_objs, xlsx_objs)

    # Save checkpoint
    saveRDS(all_objs, rds_file)
    message(">>> Saved temporary checkpoint: ", rds_file)

    # Expose to environment
    list2env(all_objs, envir = .GlobalEnv)
  }
}
