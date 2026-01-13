#' Run Correlations for Mummichog Pathway Enrichment
#'
#' Performs Pearson correlation analysis between target metabolites and the entire
#' metabolome, formatting results for use with MetaboAnalystR mummichog pathway
#' enrichment analysis.
#'
#' @param data Data frame containing metabolite intensities with ID and stage_bin columns
#' @param target_metabolites Character vector of feature names to correlate
#' @param feature_lookup Named vector mapping feature names to display names
#'
#' @return Named list of data frames, one per target metabolite, containing:
#'   - feature: Feature identifier
#'   - m.z: Mass-to-charge ratio
#'   - r.t: Retention time
#'   - p.value: Correlation p-value
#'   - statistic: Pearson correlation coefficient
#'   - method: "pearson"
#'
#' @examples
#' \dontrun{
#' results <- mummichog_correlations(
#'   data = UFT_full,
#'   target_metabolites = c("HILIC_287.236976_24.585", "C18_153.019353_22.269"),
#'   feature_lookup = c("HILIC_287.236976_24.585" = "11-cis-Retinol")
#' )
#' }
#'
#' @export
mummichog_correlations <- function(data, target_metabolites, feature_lookup) {
  # Load required libraries
  library(dplyr)
  library(purrr)
  library(stringr)
  
  # Get all metabolite columns once (exclude ID and stage_bin)
  metabolite_cols <- setdiff(names(data), c("ID", "stage_bin"))
  
  # Extract metabolite matrix for vectorized operations
  metabolite_matrix <- as.matrix(data[, metabolite_cols])
  
  # For each target metabolite, correlate it against all metabolites
  total_metabolites <- length(target_metabolites)
  results <- imap(target_metabolites, function(target_feature, idx) {
    # Progress message
    metabolite_name <- feature_lookup[target_feature]
    cat(sprintf("Processing metabolite %d of %d: %s\n", idx, total_metabolites, metabolite_name))
    
    # Extract target metabolite values
    target_values <- data[[target_feature]]
    
    # Vectorized correlation computation
    cor_coefficients <- cor(target_values, metabolite_matrix, use = "pairwise.complete.obs")
    
    # Calculate p-values using correlation test statistic
    # t = r * sqrt(n-2) / sqrt(1-r^2), df = n-2
    n <- sum(!is.na(target_values))
    t_stat <- cor_coefficients * sqrt(n - 2) / sqrt(1 - cor_coefficients^2)
    p_values <- 2 * pt(abs(t_stat), df = n - 2, lower.tail = FALSE)
    
    # Extract m/z and retention time from feature names
    mz_rt_matrix <- str_match(metabolite_cols, "(.+)_(\\d+\\.\\d+)_(\\d+\\.\\d+)")
    
    # Build results tibble
    cor_results <- tibble(
      feature = metabolite_cols,
      m.z = as.numeric(mz_rt_matrix[, 3]),
      r.t = as.numeric(mz_rt_matrix[, 4]),
      p.value = as.vector(p_values),
      statistic = as.vector(cor_coefficients),
      method = "pearson"
    )
    
    return(cor_results)
  })
  
  # Name results using feature_lookup
  names(results) <- feature_lookup[target_metabolites]
  
  return(results)
}
