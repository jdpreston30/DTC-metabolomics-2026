#' Perform t-tests between groups and format results for Mummichog analysis
#'
#' @param data Data frame with grouping column (e.g., PGD) and feature columns
#' @param group_column Name of the grouping column (default: "PGD")
#' @param output_filename Filename for the exported CSV results (optional - if NULL, no CSV is written)
#' @param output_dir Directory path for output (default: "Outputs/mummichog_inputs/")
#' @param group1_value Value representing the first group (default: "N")
#' @param group2_value Value representing the second group (default: "Y")
#' @return List containing the results tibble and summary statistics
#' @export
mummichog_ttests <- function(data,
                             group_column = "PGD",
                             output_filename = NULL,
                             output_dir = "Outputs/mummichog_inputs/",
                             group1_value = "N",
                             group2_value = "Y") {
  # _Prepare data for t-tests - data already has grouping column
  ttest_data <- data |>
    dplyr::rename(Group_Test = !!rlang::sym(group_column)) |>
    dplyr::filter(!is.na(Group_Test))

  # _Get feature names (exclude Group_Test, Patient, and any other non-numeric columns)
  excluded_columns <- c("Group_Test", "Patient")
  feature_names <- names(ttest_data)[!names(ttest_data) %in% excluded_columns]

  # _Filter to only include columns that look like metabolite features (HILIC or C18 prefix)
  metabolite_feature_names <- feature_names[stringr::str_starts(feature_names, "HILIC|C18")]

  # _Only keep numeric columns for t-tests
  numeric_feature_names <- c()
  for (col in metabolite_feature_names) {
    test_values <- ttest_data[[col]]
    if (is.numeric(test_values) || (!is.factor(test_values) && !all(is.na(suppressWarnings(as.numeric(as.character(test_values))))))) {
      numeric_feature_names <- c(numeric_feature_names, col)
    }
  }
  feature_names <- numeric_feature_names

  cat("Total columns in data:", ncol(ttest_data), "\n")
  cat("Metabolite feature columns (HILIC/C18):", length(metabolite_feature_names), "\n")
  cat("Numeric feature columns for t-tests:", length(feature_names), "\n")

  # _Initialize results list
  ttest_results <- list()

  # _Perform t-test for each feature between groups
  cat("Performing t-tests for", length(feature_names), "features...\n")

  for (i in seq_along(feature_names)) {
    feature <- feature_names[i]

    # _Progress indicator
    if (i %% 5000 == 0) {
      cat("Processed", i, "of", length(feature_names), "features...\n")
    }

    group1_values <- ttest_data[ttest_data$Group_Test == group1_value, feature]
    group2_values <- ttest_data[ttest_data$Group_Test == group2_value, feature]

    # _Convert to numeric if factor and remove NA values
    if (is.factor(group1_values)) group1_values <- as.numeric(as.character(group1_values))
    if (is.factor(group2_values)) group2_values <- as.numeric(as.character(group2_values))
    group1_values <- group1_values[!is.na(group1_values)]
    group2_values <- group2_values[!is.na(group2_values)]

    # _Check if there's sufficient data for t-test
    if (length(group1_values) < 2 || length(group2_values) < 2) {
      ttest_results[[feature]] <- 1.0 # Assign p-value of 1 for insufficient data
      next
    }

    # _Check for constant data (no variance) - assign p-value of 1
    # _Use a safer method to check for constant data
    group1_constant <- length(unique(group1_values)) == 1
    group2_constant <- length(unique(group2_values)) == 1
    if (group1_constant && group2_constant) {
      ttest_results[[feature]] <- 1.0 # Assign p-value of 1 for constant data
      next
    }

    # _Perform t-test with error handling
    tryCatch(
      {
        test_result <- t.test(group1_values, group2_values)
        ttest_results[[feature]] <- test_result$p.value
      },
      error = function(e) {
        cat("Warning: Assigning p-value of 1 to feature", feature, "due to error:", e$message, "\n")
        ttest_results[[feature]] <- 1.0 # Assign p-value of 1 for errors
      }
    )
  }

  # _Parse feature names and create results tibble
  results_tibble <- tibble::tibble(
    Feature = names(ttest_results),
    p.value = unlist(ttest_results)
  ) |>
    dplyr::mutate(
      # _Extract mode and convert: HILIC -> pos, C18 -> neg
      mode = dplyr::case_when(
        stringr::str_starts(Feature, "HILIC") ~ "positive",
        stringr::str_starts(Feature, "C18") ~ "negative",
        TRUE ~ NA_character_
      ),
      # _Extract m.z (first number after underscore)
      m.z = stringr::str_extract(Feature, "(?<=_)[0-9.]+"),
      # _Extract rt (second number - after second underscore)
      r.t = stringr::str_extract(Feature, "_[0-9.]+_([0-9.]+)") |>
        stringr::str_extract("[0-9.]+$")
    ) |>
    dplyr::select(m.z, p.value, mode, r.t) |>
    dplyr::mutate(
      m.z = as.numeric(m.z),
      r.t = as.numeric(r.t)
    ) |>
    # _Remove rows where feature parsing failed (invalid feature names)
    dplyr::filter(!is.na(m.z) & !is.na(mode))

  # _Export results (only if output_filename is provided)
  output_path <- NULL
  if (!is.null(output_filename)) {
    # _Create output directory if it doesn't exist
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
      cat("Created output directory:", output_dir, "\n")
    }

    # _Create full output path
    output_path <- file.path(output_dir, output_filename)

    # _Export results
    readr::write_csv(results_tibble, output_path)
    cat("T-test results exported to:", output_path, "\n")
  } else {
    cat("No output filename provided - results returned without writing CSV\n")
  }

  # _Display summary
  cat("T-test results exported to:", output_path, "\n")
  cat("Total features processed:", length(feature_names), "\n")
  cat("Features with actual t-test p-values:", sum(unlist(ttest_results) < 1.0), "\n")
  cat("Features with assigned p-value of 1 (constant/insufficient data):", sum(unlist(ttest_results) == 1.0), "\n")
  cat("Total results:", nrow(results_tibble), "\n")
  cat("Group", group1_value, "vs Group", group2_value, "comparison\n")
  cat("First 10 results:\n")
  print(head(results_tibble, 10))

  # _Return results
  return(list(
    results = results_tibble,
    n_features = nrow(results_tibble),
    output_file = output_path, # Will be NULL if no CSV written
    group1_value = group1_value,
    group2_value = group2_value
  ))
}