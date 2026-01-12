#' Remove Low-Variance Features Based on Identical Value Percentage
#'
#' This function removes feature columns that have a percentage of identical values
#' greater than a specified threshold, which may indicate poor data quality or
#' QC issues. Metadata columns (ID, stage_bin) are automatically excluded from
#' the filtering and retained in the output.
#'
#' @param data A data frame or tibble containing feature data
#' @param threshold Numeric value between 0 and 1 representing the maximum allowed
#'   percentage of identical values. Columns with a higher percentage will be removed.
#'   Default is 0.2 (20%).
#' @param metadata_cols Character vector of metadata column names to exclude from
#'   filtering (default c("ID", "stage_bin"))
#'
#' @return A tibble with low-variance feature columns removed
#'
#' @examples
#' \dontrun{
#' # Remove features with >20% identical values (default)
#' clean_data <- remove_qc(UFT_filtered)
#'
#' # Use stricter threshold: remove features with >10% identical values
#' clean_data <- remove_qc(UFT_filtered, threshold = 0.1)
#'
#' # More lenient: remove only features with >50% identical values
#' clean_data <- remove_qc(UFT_filtered, threshold = 0.5)
#' }
#'
#' @export
remove_qc <- function(data, threshold = 0.2, metadata_cols = c("ID", "stage_bin")) {
  library(dplyr)
  library(tidyr)
  
  # Identify feature columns (exclude metadata)
  feature_cols <- setdiff(names(data), metadata_cols)
  
  # Calculate percentage of identical values for each feature column
  identical_pcts <- data |>
    select(all_of(feature_cols)) |>
    summarise(across(everything(), ~max(table(.)) / length(.))) |>
    pivot_longer(everything(), names_to = "feature", values_to = "max_identical_pct")
  
  # Identify columns to keep (below or equal to threshold)
  cols_to_keep <- identical_pcts |>
    filter(max_identical_pct <= threshold) |>
    pull(feature)
  
  # Count removed columns
  n_removed <- length(feature_cols) - length(cols_to_keep)
  
  # Report
  cat(
    "\n", strrep("=", 60), "\n",
    "QC FILTER: Removing Low-Variance Features\n",
    strrep("=", 60), "\n\n",
    "Threshold: ", threshold * 100, "% identical values\n",
    "Features analyzed: ", length(feature_cols), "\n",
    "Features removed: ", n_removed, "\n",
    "Features retained: ", length(cols_to_keep), "\n\n",
    strrep("=", 60), "\n\n"
  )
  
  # Return data with only retained columns plus metadata
  data |>
    select(all_of(c(metadata_cols, cols_to_keep)))
}
