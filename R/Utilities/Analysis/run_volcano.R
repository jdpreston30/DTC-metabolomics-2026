#' Run volcano plot analysis with automatic t-tests on feature table
#'
#' @param data Data frame with patient IDs, grouping variable, and feature columns
#' @param group_var Character string specifying the column name to use for grouping
#' @param patient_var Character string specifying the Patient ID column name (default: "Patient")
#' @param group_levels Optional vector specifying factor levels for grouping (must be exactly 2 levels)
#' @param fc_threshold Fold change threshold for significance (default: log2(1.5) = 0.585)
#' @param p_threshold P-value threshold for significance (default: 0.05)
#' @return List containing volcano_data and analysis metadata
#' @export
run_volcano <- function(data,
                        group_var,
                        patient_var = "Patient",
                        group_levels = NULL,
                        fc_threshold = log2(1.5),
                        p_threshold = 0.05) {
  # ---- Data preparation ----
  if (!group_var %in% names(data)) {
    stop(paste("Group variable", group_var, "not found in data"))
  }
  if (!patient_var %in% names(data)) {
    stop(paste("Patient variable", patient_var, "not found in data"))
  }

  dat <- as.data.frame(data)

  # Set up grouping variable
  if (!is.null(group_levels)) {
    dat[[group_var]] <- factor(dat[[group_var]], levels = group_levels)
  } else {
    dat[[group_var]] <- factor(dat[[group_var]])
  }

  # Check for exactly 2 groups
  n_groups <- length(levels(dat[[group_var]]))
  if (n_groups != 2) {
    stop(paste("Volcano plot requires exactly 2 groups, but found", n_groups, "groups"))
  }

  # Identify numeric columns (exclude patient and grouping variables)
  factor_cols <- sapply(dat, is.factor)
  cols_to_exclude <- c(patient_var, group_var)
  other_factor_cols <- names(factor_cols)[factor_cols & !names(factor_cols) %in% cols_to_exclude]
  drop_cols <- c(cols_to_exclude, other_factor_cols)
  numeric_cols <- !names(dat) %in% drop_cols

  if (sum(numeric_cols) == 0) {
    stop("No numeric columns found for t-tests")
  }

  feature_data <- dat[, numeric_cols, drop = FALSE]
  group <- dat[[group_var]]

  message(
    "Performing t-tests on ", ncol(feature_data), " features between ",
    paste(levels(group), collapse = " vs ")
  )

  # ---- Perform t-tests ----
  ttest_results <- data.frame(
    feature = colnames(feature_data),
    p_value = NA,
    log2_fc = NA,
    mean_group1 = NA,
    mean_group2 = NA,
    stringsAsFactors = FALSE
  )

  group_names <- levels(group)

  for (i in seq_len(ncol(feature_data))) {
    feature_values <- feature_data[, i]

    # Skip if all values are identical
    if (length(unique(feature_values)) == 1) {
      next
    }

    # Calculate group means (log2 scale)
    group1_vals <- feature_values[group == group_names[1]]
    group2_vals <- feature_values[group == group_names[2]]

    mean1 <- mean(group1_vals, na.rm = TRUE)
    mean2 <- mean(group2_vals, na.rm = TRUE)

    # Calculate log2 fold change (group2 vs group1)
    log2_fc <- mean2 - mean1

    # Perform t-test
    tryCatch(
      {
        t_result <- t.test(feature_values ~ group)
        ttest_results$p_value[i] <- t_result$p.value
        ttest_results$log2_fc[i] <- log2_fc
        ttest_results$mean_group1[i] <- mean1
        ttest_results$mean_group2[i] <- mean2
      },
      error = function(e) {
        # Skip features that cause t-test errors
      }
    )
  }

  # Remove failed tests
  ttest_results <- ttest_results[!is.na(ttest_results$p_value), ]

  # ---- Create volcano data ----
  volcano_data <- ttest_results
  volcano_data$neg_log10_p <- -log10(volcano_data$p_value)

  # Classify significance with custom legend text
  volcano_data$Legend <- "Not Significant"

  # Determine legend labels based on comparison type
  if ("Early" %in% group_names && "Advanced" %in% group_names) {
    # Mild/Mod. PGD vs Severe PGD comparison: both are "Up in" labels
    up_label <- "Up in Stage III/IV"
    down_label <- "Down in Stage III/IV"
  } else {
    # Other comparisons (e.g., No PGD vs Severe): traditional up/down
    up_label <- "Up in Stage III/IV"
    down_label <- "Down in Stage III/IV"
  }

  volcano_data$Legend[volcano_data$p_value < p_threshold & volcano_data$log2_fc > fc_threshold] <- up_label
  volcano_data$Legend[volcano_data$p_value < p_threshold & volcano_data$log2_fc < -fc_threshold] <- down_label

  volcano_data$Legend <- factor(volcano_data$Legend,
    levels = c("Not Significant", up_label, down_label)
  )

  # Print summary statistics
  cat("\n=== Volcano Analysis Summary ===\n")
  cat("Highest -log10(p-value):", round(max(volcano_data$neg_log10_p, na.rm = TRUE), 3), "\n")
  cat("Min log2(FC):", round(min(volcano_data$log2_fc, na.rm = TRUE), 3), "\n")
  cat("Max log2(FC):", round(max(volcano_data$log2_fc, na.rm = TRUE), 3), "\n")
  cat("================================\n\n")

  # Return analysis results
  list(
    volcano_data = volcano_data,
    group_names = group_names,
    group_var = group_var,
    patient_var = patient_var,
    fc_threshold = fc_threshold,
    p_threshold = p_threshold,
    up_label = up_label,
    down_label = down_label,
    original_data = dat
  )
}
