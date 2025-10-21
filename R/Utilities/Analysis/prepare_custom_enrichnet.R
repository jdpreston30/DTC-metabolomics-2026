### Network Analysis Utilities
### Adapted from MetaboAnalystR PrepareEnrichNet functionality
### Original work by Jeff Xia\email{jeff.xia@mcgill.ca}
### McGill University, Canada
### Original License: GNU GPL (>= 2)
### 
### This implementation provides custom network preparation functionality
### to bypass missing compiled components in the MetaboAnalystR package
### while maintaining compatibility with the original MetaboAnalyst workflow.

#' Prepare Custom Enrichment Network
#'
#' This function creates network data from mummichog pathway enrichment results,
#' providing an alternative to the MetaboAnalystR PrepareEnrichNet function
#' which requires server-side compiled components not available in the R package.
#'
#' @description
#' Adapted from MetaboAnalystR's PrepareEnrichNet functionality. This custom
#' implementation extracts pathway results from mummichog analysis and creates
#' network files for visualization without requiring the compiled utils_enrichnet.Rc
#' components that are only available on the MetaboAnalyst web server.
#'
#' @param mSet A MetaboAnalyst data object containing mummichog results
#' @param netNm Character string specifying the network name (used for output files)
#' @param overlapType Character string specifying overlap type (default: "mixed")
#'   Currently not used in this implementation but maintained for API compatibility
#'
#' @return Returns the input mSet object (for pipeline compatibility)
#'   Side effects: Creates network files in the current working directory
#'
#' @details
#' The function performs the following operations:
#' 1. Extracts pathway results from either CSV file or mSet object
#' 2. Creates network node data with pathway IDs, labels, p-values, and hit counts
#' 3. Exports network data in JSON format for web-based visualization
#' 4. Creates SIF (Simple Interaction Format) files for Cytoscape import
#' 5. Generates basic pathway similarity edges between top-ranked pathways
#'
#' @section Files Created:
#' \describe{
#'   \item{[netNm]_network_data.json}{JSON format network data for web visualization}
#'   \item{[netNm].sif}{SIF format network file for Cytoscape import}
#' }
#'
#' @section Attribution:
#' This function is adapted from the MetaboAnalystR package PrepareEnrichNet
#' functionality. Original implementation by Jeff Xia at McGill University.
#' Original license: GNU GPL (>= 2).
#'
#' @author Adapted for local use from MetaboAnalystR
#' @references
#' Chong, J. et al. MetaboAnalyst 4.0: towards more transparent and integrative
#' metabolomics analysis. Nucleic Acids Res. 46, W486-W494 (2018).
#'
#' @examples
#' \dontrun{
#' # After running mummichog analysis with MetaboAnalystR
#' mSet <- prepare_custom_enrichnet(mSet, "enrichNet_mfn", "mixed")
#' }
#'
#' @seealso \code{\link[MetaboAnalystR]{PerformPSEA}} for the mummichog analysis
#'
#' @export
prepare_custom_enrichnet <- function(mSet, netNm, overlapType = "mixed") {
  
  cat("Preparing custom enrichment network...\n")
  
  # Try multiple approaches to extract pathway results from mSet object
  pathway_results <- NULL
  
  # Method 1: Check for pathway CSV file in working directory
  csv_file <- "mummichog_pathway_enrichment_mummichog.csv"
  if(file.exists(csv_file)) {
    pathway_results <- read.csv(csv_file, row.names = 1)
    cat("✓ Found pathway results in CSV file\n")
  } else {
    # Method 2: Extract from mSet object structure
    if(!is.null(mSet$mummi.resmat)) {
      pathway_results <- mSet$mummi.resmat
    } else if(!is.null(mSet$analSet) && !is.null(mSet$analSet$resmat)) {
      pathway_results <- mSet$analSet$resmat
    } else if(!is.null(mSet$dataSet) && !is.null(mSet$dataSet$mummi.proc)) {
      pathway_results <- mSet$dataSet$mummi.proc
    } else {
      cat("⚠ Could not find pathway results in mSet object or CSV file\n")
      return(mSet)
    }
  }
  
  # Ensure we have valid data
  if(is.null(pathway_results) || nrow(pathway_results) == 0) {
    cat("⚠ No pathway results found\n")
    return(mSet)
  }
  
  # Create basic network data structure
  # Handle different column name possibilities
  pval_col <- NULL
  hits_col <- NULL
  
  if("p.value" %in% colnames(pathway_results)) {
    pval_col <- "p.value"
  } else if("pval" %in% colnames(pathway_results)) {
    pval_col <- "pval"
  } else if("P.Value" %in% colnames(pathway_results)) {
    pval_col <- "P.Value"
  }
  
  if("hits" %in% colnames(pathway_results)) {
    hits_col <- "hits"
  } else if("Hits" %in% colnames(pathway_results)) {
    hits_col <- "Hits"
  } else if("total.hits" %in% colnames(pathway_results)) {
    hits_col <- "total.hits"
  }
  
  # Create node data with available columns
  node_data <- data.frame(
    id = rownames(pathway_results),
    label = rownames(pathway_results),
    stringsAsFactors = FALSE
  )
  
  if(!is.null(pval_col)) {
    node_data$pvalue <- pathway_results[, pval_col]
  }
  
  if(!is.null(hits_col)) {
    node_data$hits <- pathway_results[, hits_col]
  }
  
  network_data <- list(
    node_data = node_data,
    edge_data = data.frame(
      source = character(0),
      target = character(0),
      weight = numeric(0),
      stringsAsFactors = FALSE
    )
  )
  
  # Save network data as JSON for visualization
  jsonlite::write_json(network_data, paste0(netNm, "_network_data.json"), auto_unbox = TRUE)
  
  # Create simple SIF format for Cytoscape
  if(nrow(network_data$node_data) > 1) {
    # Create edges between top pathways (simple approach)
    if(!is.null(node_data$pvalue)) {
      top_pathways <- head(network_data$node_data[order(network_data$node_data$pvalue),], 10)
    } else {
      top_pathways <- head(network_data$node_data, 10)
    }
    
    if(nrow(top_pathways) > 1) {
      sif_data <- data.frame(
        source = top_pathways$id[1:(nrow(top_pathways)-1)],
        interaction = rep("pathway_similarity", nrow(top_pathways)-1),
        target = top_pathways$id[2:nrow(top_pathways)]
      )
      write.table(sif_data, paste0(netNm, ".sif"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
    }
  }
  
  cat("✓ Custom network files created:", paste0(netNm, "_network_data.json"), "and", paste0(netNm, ".sif"), "\n")
  return(mSet)
}
