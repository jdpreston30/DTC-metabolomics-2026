#' Run T-tests on Targeted Metabolomic Features
#'
#' This function performs t-tests for all C18 and HILIC features in a feature table
#' against a grouping variable, and returns comprehensive results including means,
#' p-values, FDR correction, and feature annotations.
#'
#' @param feature_table Data frame containing feature data with samples as rows and features as columns
#' @param grouping_var Either a vector indicating group membership OR a column name in feature_table
#' @param p_adjust_method Method for p-value adjustment (default "fdr")
#' @param fc_ref_group Reference group for fold change calculation (optional). If specified, fold change will be calculated as other_group/ref_group using raw intensities (2^log2_mean)
#'
#' @return A tibble with columns:
#'   - feature: Feature identifier
#'   - display_name: Feature name from TFT_annot key
#'   - mean_[group1]: Mean in first factor level
#'   - mean_[group2]: Mean in second factor level
#'   - mean_overall: Overall population mean (combining both groups)
#'   - fold_change: Fold change calculated as other_group/ref_group using raw intensities (if fc_ref_group specified)
#'   - p_value: Raw p-value from t-test
#'   - p_value_fdr: FDR-corrected p-value
#'   - unique_vals_no_severe: Unique values in FALSE group (n_unique/n_total)
#'   - unique_vals_severe: Unique values in TRUE group (n_unique/n_total)
#'   - unique_vals: Overall unique values (n_unique/n_total)
#'   - unique_percentage: Percentage of unique values in the dataset
#'   - low_detect_likely: "Y" if unique_percentage > 20%, "N" otherwise
#'   - isomer: Isomer status from TFT_annot key
#'   - MMD: Multi-mode detection status from TFT_annot key
#'   - n_[group1]: Sample size in first group
#'   - n_[group2]: Sample size in second group
#'
#' @examples
#' \dontrun{
#'   # Using column name
#'   results <- run_targeted_ttests(
#'     feature_table = TFT_annot,
#'     grouping_var = "severe_PGD",
#'     fc_ref_group = "No Severe PGD"
#'   )
#'   
#'   # Using vector
#'   results <- run_targeted_ttests(
#'     feature_table = TFT_annot,
#'     grouping_var = clinical_metadata$Severe_PGD_factor
#'   )
#' }
#'
#' @export
run_targeted_ttests <- function(feature_table, 
                               grouping_var,
                               p_adjust_method = "fdr",
                               fc_ref_group = NULL) {
  
  # Load required libraries
  library(dplyr)
  library(purrr)
  library(broom)
  
  # Handle grouping_var - could be column name or vector
  if (is.character(grouping_var) && length(grouping_var) == 1) {
    # It's a column name
    if (!grouping_var %in% colnames(feature_table)) {
      stop("Column '", grouping_var, "' not found in feature_table")
    }
    group_vector <- feature_table[[grouping_var]]
  } else {
    # It's a vector
    group_vector <- grouping_var
  }
  
  # Convert to factor if not already
  if (!is.factor(group_vector)) {
    group_vector <- as.factor(group_vector)
    cat("Converted grouping variable to factor\n")
  }
  
  # Validate factor has 2 levels
  if (length(levels(group_vector)) != 2) {
    stop("grouping_var must have exactly 2 levels, found: ", 
         paste(levels(group_vector), collapse = ", "))
  }
  
  if (nrow(feature_table) != length(group_vector)) {
    stop("Number of rows in feature_table (", nrow(feature_table), 
         ") must match length of grouping_var (", length(group_vector), ")")
  }
  
  # Get feature names that start with C18 or HILIC
  feature_names <- colnames(feature_table)
  metabolomic_features <- feature_names[grepl("^(C18_|HILIC_)", feature_names)]
  
  if (length(metabolomic_features) == 0) {
    stop("No features starting with 'C18_' or 'HILIC_' found in feature_table")
  }
  
  cat("Found", length(metabolomic_features), "metabolomic features to test\n")
  
  # Get group levels
  group_levels <- levels(group_vector)
  group1_name <- group_levels[1]
  group2_name <- group_levels[2]
  
  cat("Comparing:", group1_name, "vs", group2_name, "\n")
  
  # Function to perform t-test for a single feature
  perform_ttest <- function(feature_name) {
    # Get feature values
    feature_values <- feature_table[[feature_name]]
    
    # Remove any missing values
    complete_cases <- !is.na(feature_values) & !is.na(group_vector)
    feature_clean <- feature_values[complete_cases]
    group_clean <- group_vector[complete_cases]
    
    # Check if we have enough data
    if (length(feature_clean) < 3 || length(unique(group_clean)) < 2) {
      tibble(
        feature = feature_name,
        mean_group1 = NA_real_,
        mean_group2 = NA_real_,
        mean_overall = NA_real_,
        fold_change = NA_real_,
        p_value = NA_real_,
        unique_vals_no_severe = "0/0",
        unique_vals_severe = "0/0", 
        unique_vals = "0/0",
        unique_percentage = 0.0,
        low_detect_likely = "N",
        n_group1 = 0,
        n_group2 = 0
      )
    }
    
    # Calculate group means
    group1_values <- feature_clean[group_clean == group1_name]
    group2_values <- feature_clean[group_clean == group2_name]
    
    mean_group1 <- mean(group1_values, na.rm = TRUE)
    mean_group2 <- mean(group2_values, na.rm = TRUE)
    mean_overall <- mean(feature_clean, na.rm = TRUE)
    
    # Calculate fold change if fc_ref_group is specified
    fold_change <- NA_real_
    if (!is.null(fc_ref_group)) {
      # Convert log2 transformed means back to raw intensities
      raw_mean_group1 <- 2^mean_group1
      raw_mean_group2 <- 2^mean_group2
      
      # Calculate fold change based on reference group
      if (fc_ref_group == group1_name) {
        # group2 / group1 (fold change relative to group1)
        fold_change <- raw_mean_group2 / raw_mean_group1
      } else if (fc_ref_group == group2_name) {
        # group1 / group2 (fold change relative to group2)
        fold_change <- raw_mean_group1 / raw_mean_group2
      } else {
        warning("fc_ref_group '", fc_ref_group, "' not found in group levels. Available: ", 
                paste(group_levels, collapse = ", "))
      }
    }
    
    # Calculate unique values analysis
    # Get all values for this feature across the entire dataset for duplicate detection
    all_feature_values <- feature_table[[feature_name]][!is.na(feature_table[[feature_name]])]
    
    # Function to count unique values (values that appear only once in the entire dataset)
    count_unique_in_subgroup <- function(subgroup_values, all_values) {
      # Count how many times each value appears in the entire dataset
      value_counts <- table(all_values)
      # In the subgroup, count how many values appear only once in the entire dataset
      unique_in_subgroup <- sum(sapply(subgroup_values, function(x) value_counts[as.character(x)] == 1))
      return(unique_in_subgroup)
    }
    
    # Calculate unique values for each subgroup
    unique_count_group1 <- count_unique_in_subgroup(group1_values, all_feature_values)
    unique_count_group2 <- count_unique_in_subgroup(group2_values, all_feature_values)
    unique_count_total <- count_unique_in_subgroup(feature_clean, all_feature_values)
    
    # Create formatted strings (assuming group1 is FALSE/no severe, group2 is TRUE/severe)
    unique_vals_no_severe <- paste0(unique_count_group1, "/", length(group1_values))
    unique_vals_severe <- paste0(unique_count_group2, "/", length(group2_values))
    unique_vals <- paste0(unique_count_total, "/", length(feature_clean))
    
    # Calculate percentage of unique values
    unique_percentage <- (unique_count_total / length(feature_clean)) * 100

    # Determine if low detection is likely (>80% unique values)
    low_detect_likely <- if (unique_percentage > 80) "N" else "Y"
    
    # Perform t-test
    tryCatch({
      t_result <- t.test(feature_clean ~ group_clean)
      
      tibble(
        feature = feature_name,
        mean_group1 = mean_group1,
        mean_group2 = mean_group2,
        mean_overall = mean_overall,
        fold_change = fold_change,
        log2FC = log2(fold_change),
        p_value = t_result$p.value,
        unique_vals_no_severe = unique_vals_no_severe,
        unique_vals_severe = unique_vals_severe,
        unique_vals = unique_vals,
        unique_percentage = round(unique_percentage, 1),
        low_detect_likely = low_detect_likely,
        n_group1 = length(group1_values),
        n_group2 = length(group2_values)
      )
    }, error = function(e) {
      tibble(
        feature = feature_name,
        mean_group1 = mean_group1,
        mean_group2 = mean_group2,
        mean_overall = mean_overall,
        fold_change = fold_change,
        log2FC = log2(fold_change),
        p_value = NA_real_,
        unique_vals_no_severe = unique_vals_no_severe,
        unique_vals_severe = unique_vals_severe,
        unique_vals = unique_vals,
        unique_percentage = round(unique_percentage, 1),
        low_detect_likely = low_detect_likely,
        n_group1 = length(group1_values),
        n_group2 = length(group2_values)
      )
    })
  }
  
  # Run t-tests for all features
  cat("Running t-tests...\n")
  ttest_results <- map_dfr(metabolomic_features, perform_ttest)
  
  # Apply FDR correction
  ttest_results <- ttest_results %>%
    mutate(
      p_value_fdr = p.adjust(p_value, method = p_adjust_method)
    ) %>%
    # Reorder columns to put p_value_fdr right after p_value
    select(feature:p_value, p_value_fdr, everything())
  
  # Join with TFT_annot key to get feature annotations (columns already in snake_case)
  final_results <- ttest_results %>%
    arrange(p_value)
  
  # Convert group names to snake_case for column names
  group1_snake <- tolower(gsub("[^A-Za-z0-9]", "_", group1_name))
  group2_snake <- tolower(gsub("[^A-Za-z0-9]", "_", group2_name))
  
  # Add group names to column names for clarity (snake_case)
  colnames(final_results)[colnames(final_results) == "mean_group1"] <- paste0("mean_", group1_snake)
  colnames(final_results)[colnames(final_results) == "mean_group2"] <- paste0("mean_", group2_snake)
  colnames(final_results)[colnames(final_results) == "n_group1"] <- paste0("n_", group1_snake)
  colnames(final_results)[colnames(final_results) == "n_group2"] <- paste0("n_", group2_snake)
  
  # Summary statistics
  n_tested <- nrow(final_results)
  n_significant_raw <- sum(final_results$p_value < 0.05, na.rm = TRUE)
  n_significant_fdr <- sum(final_results$p_value_fdr < 0.05, na.rm = TRUE)
  
  cat("\nResults Summary:\n")
  cat("- Total features tested:", n_tested, "\n")
  cat("- Significant at p < 0.05:", n_significant_raw, "\n")
  cat("- Significant at FDR < 0.05:", n_significant_fdr, "\n")
  
  return(final_results)
}
