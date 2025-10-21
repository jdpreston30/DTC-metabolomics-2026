#' Create Biological Network Based on Shared Metabolites
#'
#' This function creates a biologically meaningful network by connecting pathways
#' that share empirical compounds (metabolites) detected in the mummichog analysis.
#' This recreates what the real PrepareEnrichNet function does on the MetaboAnalyst server.
#'
#' @param pathway_csv Path to the mummichog pathway enrichment CSV file
#' @param min_shared_compounds Minimum number of shared compounds to create an edge (default: 2)
#' @param p_threshold P-value threshold for including pathways (default: 0.1)
#' @param max_pathways Maximum number of pathways to include (default: 20)
#' @param network_name Name for output files (default: "biological_network")
#'
#' @return List containing network data and creates SIF and JSON files
#'
#' @export
create_biological_network <- function(pathway_csv, min_shared_compounds = 2, 
                                     p_threshold = 0.1, max_pathways = 20,
                                     network_name = "biological_network") {
  
  library(dplyr)
  library(jsonlite)
  
  cat("Creating biological network from pathway results...\n")
  
  # Read pathway results
  if (!file.exists(pathway_csv)) {
    stop("Pathway CSV file not found: ", pathway_csv)
  }
  
  pathway_data <- read.csv(pathway_csv, stringsAsFactors = FALSE, row.names = 1)
  
  # Filter significant pathways and limit to top pathways
  significant_pathways <- pathway_data %>%
    filter(P.Fisher. <= p_threshold) %>%
    arrange(P.Fisher.) %>%
    head(max_pathways)
  
  if (nrow(significant_pathways) == 0) {
    stop("No significant pathways found with p <= ", p_threshold)
  }
  
  cat("Found", nrow(significant_pathways), "significant pathways\n")
  
  # Debug: Check what pathway names we have
  cat("Pathway names:", paste(head(rownames(significant_pathways)), collapse = ", "), "\n")
  
  # Parse compound hits for each pathway
  pathway_compounds <- list()
  for (i in 1:nrow(significant_pathways)) {
    pathway_name <- rownames(significant_pathways)[i]
    compounds_str <- significant_pathways$cpd.hits[i]
    
    # Split compound string by semicolon
    if (!is.na(compounds_str) && compounds_str != "") {
      compounds <- trimws(strsplit(compounds_str, ";")[[1]])
      pathway_compounds[[pathway_name]] <- compounds
    } else {
      pathway_compounds[[pathway_name]] <- character(0)
    }
  }
  
  # Calculate shared compounds between all pathway pairs
  pathway_names <- names(pathway_compounds)
  edges_list <- list()
  
  for (i in 1:(length(pathway_names)-1)) {
    for (j in (i+1):length(pathway_names)) {
      pathway1 <- pathway_names[i]
      pathway2 <- pathway_names[j]
      
      compounds1 <- pathway_compounds[[pathway1]]
      compounds2 <- pathway_compounds[[pathway2]]
      
      # Find shared compounds
      shared_compounds <- intersect(compounds1, compounds2)
      n_shared <- length(shared_compounds)
      
      # Create edge if enough shared compounds
      if (n_shared >= min_shared_compounds) {
        edges_list[[length(edges_list) + 1]] <- data.frame(
          source = pathway1,
          target = pathway2,
          interaction = "shared_metabolites",
          weight = n_shared,
          shared_compounds = paste(shared_compounds, collapse = ";"),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  
  # Combine all edges
  if (length(edges_list) > 0) {
    edges_df <- do.call(rbind, edges_list)
    cat("Created", nrow(edges_df), "biological edges based on shared metabolites\n")
  } else {
    edges_df <- data.frame(
      source = character(0),
      target = character(0), 
      interaction = character(0),
      weight = numeric(0),
      shared_compounds = character(0),
      stringsAsFactors = FALSE
    )
    cat("No biological edges found (try lowering min_shared_compounds)\n")
  }
  
  # Create node data
  nodes_df <- data.frame(
    id = rownames(significant_pathways),
    label = rownames(significant_pathways),
    pvalue = significant_pathways$P.Fisher.,
    hits_sig = significant_pathways$Hits.sig,
    hits_total = significant_pathways$Hits.total,
    pathway_total = significant_pathways$Pathway.total,
    enrichment_factor = significant_pathways$Hits.sig / pmax(significant_pathways$Expected, 0.1),
    stringsAsFactors = FALSE
  )
  
  # Create network data structure
  network_data <- list(
    node_data = nodes_df,
    edge_data = edges_df,
    metadata = list(
      min_shared_compounds = min_shared_compounds,
      p_threshold = p_threshold,
      n_pathways = nrow(nodes_df),
      n_edges = nrow(edges_df)
    )
  )
  
  # Get the directory of the input CSV file
  csv_dir <- dirname(pathway_csv)
  
  # Save as JSON in the same directory as the CSV
  json_file <- file.path(csv_dir, paste0(network_name, "_biological_network.json"))
  write_json(network_data, json_file, auto_unbox = TRUE, pretty = TRUE)
  cat("✓ Created JSON file:", json_file, "\n")
  
  # Save as SIF (only if there are edges) in the same directory as the CSV
  if (nrow(edges_df) > 0) {
    sif_file <- file.path(csv_dir, paste0(network_name, "_biological.sif"))
    sif_data <- data.frame(
      source = edges_df$source,
      interaction = paste0("shared_metabolites_", edges_df$weight),
      target = edges_df$target
    )
    write.table(sif_data, sif_file, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
    cat("✓ Created SIF file:", sif_file, "\n")
  }
  
  # Print summary
  cat("\n=== BIOLOGICAL NETWORK SUMMARY ===\n")
  cat("Pathways included:", nrow(nodes_df), "\n")
  cat("Biological edges:", nrow(edges_df), "\n")
  cat("Based on shared metabolites with minimum", min_shared_compounds, "compounds\n")
  
  if (nrow(edges_df) > 0) {
    cat("Edge weights (shared compounds):", paste(range(edges_df$weight), collapse = " to "), "\n")
  }
  
  return(network_data)
}
