promote_checkpoint <- function(raw_path) {
  checkpoint_dir <- here::here("R", "Checkpoints")
  temp_rds <- file.path(checkpoint_dir, "all_objects.rds")
  final_rds <- file.path(raw_path, "all_objects.rds")

  if (file.exists(temp_rds)) {
    file.copy(temp_rds, final_rds, overwrite = TRUE)
    message(">>> Final checkpoint saved to: ", final_rds)
    unlink(checkpoint_dir, recursive = TRUE, force = TRUE)
    message(">>> Temporary checkpoint directory deleted: ", checkpoint_dir)
  } else {
    message("⚠️ No temporary checkpoint found at end of pipeline")
  }
}