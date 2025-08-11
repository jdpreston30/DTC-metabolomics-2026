#' Create and display heatmap with ANOVA feature selection and extract clades
#'
#' @param data Data frame with Patient_ID, Variant, and feature columns
#' @param variant_colors Named vector of colors for each variant
#' @param top_features Number of top features to include (default: 2500)
#' @param variant_levels Character vector specifying order of variants
#' @param n_clades Number of clades to extract from hierarchical clustering (default: 2)
#' @param output_filename Optional filename for saving heatmap as SVG (without path)
#' @return List containing the filtered matrix M, sample clustering results, clade assignments, and other objects
#' @export
make_heatmap <- function(data, 
                        variant_colors, 
                        top_features = 2500,
                        variant_levels = c("PTC", "FV-PTC", "FTC"),
                        n_clades = 2,
                        output_filename = NULL) {
  
  #_Break into groups
  dat <- data %>% dplyr::select(-Patient_ID)
  group <- factor(dat$Variant, levels = variant_levels)
  X <- as.matrix(dat %>% dplyr::select(-Variant)) # rows = samples, cols = features
  sample_ids <- make.names(data$Patient_ID, unique = TRUE)
  
  #_ANOVA p-values per feature (over group)
  pvals <- apply(X, 2, function(x) {
    fit <- aov(x ~ group)
    summary(fit)[[1]][["Pr(>F)"]][1]
  })
  
  #_Top features (or all if fewer)
  k <- min(top_features, ncol(X))
  top_idx <- order(pvals)[1:k]
  X_top <- X[, top_idx, drop = FALSE]
  
  #_Drop zero-variance features
  nzv <- apply(X_top, 2, sd, na.rm = TRUE) > 0
  X_top <- X_top[, nzv, drop = FALSE]
  
  #_Assign row and column names
  rownames(X_top) <- sample_ids # samples
  colnames(X_top) <- colnames(data)[-(1:2)][top_idx][nzv] # features
  M <- t(X_top)
  
  #_Column annotation (must match columns of M)
  ann_col <- data.frame(Variant = group)
  rownames(ann_col) <- sample_ids
  ann_col <- ann_col[colnames(M), , drop = FALSE]
  
  #_Custom colors for Variant
  ann_colors <- list(Variant = variant_colors)
  
  #_Create and display heatmap
  heatmap_plot <- pheatmap::pheatmap(
    M,
    scale = "row",
    color = colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(255),
    clustering_distance_rows = "euclidean",
    clustering_distance_cols = "euclidean",
    clustering_method = "complete",
    annotation_col = ann_col,
    annotation_colors = ann_colors,
    show_rownames = FALSE,
    show_colnames = FALSE,
    fontsize = 10,
    na_col = "#DDDDDD"
  )
  
  #_Save heatmap as SVG if filename provided
  if (!is.null(output_filename)) {
    svg(filename = paste0("Outputs/", output_filename), 
        width = 8, height = 8, pointsize = 12)
    pheatmap::pheatmap(
      M,
      scale = "row",
      color = colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(255),
      clustering_distance_rows = "euclidean",
      clustering_distance_cols = "euclidean",
      clustering_method = "complete",
      annotation_col = ann_col,
      annotation_colors = ann_colors,
      show_rownames = FALSE,
      show_colnames = FALSE,
      fontsize = 10,
      na_col = "#DDDDDD"
    )
    dev.off()
    cat("Heatmap saved to:", paste0("Outputs/", output_filename), "\n")
  }
  
  #_Perform sample clustering (same as pheatmap)
  Mz <- t(scale(t(M), center = TRUE, scale = TRUE))
  Mz[is.na(Mz)] <- 0
  d_cols <- dist(t(Mz), method = "euclidean")
  hc_cols <- hclust(d_cols, method = "complete")
  
  #_Extract clades from hierarchical clustering
  #_Create ID mapping for original Patient_IDs
  id_map <- setNames(data$Patient_ID, sample_ids)
  #_Cut dendrogram into clades
  clades <- cutree(hc_cols, k = n_clades)
  #_Order Patient_IDs along dendrogram
  ids_ordered_clean <- names(clades)[hc_cols$order]
  ids_ordered_orig <- unname(id_map[ids_ordered_clean])
  #_Create clade data frame
  clade_df <- tibble::tibble(
    Patient_ID = ids_ordered_orig,
    Clade = unname(clades[ids_ordered_clean])
  )
  
  #_Extract Patient_IDs for each clade
  clade_lists <- list()
  for (i in 1:n_clades) {
    clade_lists[[paste0("clade", i, "_ids")]] <- clade_df %>%
      dplyr::filter(Clade == i) %>%
      dplyr::pull(Patient_ID)
  }
  
  #_Display clade results
  cat("Patient clade assignments:\n")
  print(clade_df)
  for (i in 1:n_clades) {
    cat("Clade", i, "patients:", length(clade_lists[[paste0("clade", i, "_ids")]]), "\n")
  }
  
  #_Return useful objects
  return(list(
    M = M,
    Mz = Mz,
    hc_cols = hc_cols,
    sample_ids = sample_ids,
    group = group,
    pvals = pvals,
    top_idx = top_idx,
    nzv = nzv,
    heatmap_plot = heatmap_plot,
    clade_df = clade_df,
    clades = clades,
    clade_lists = clade_lists,
    id_map = id_map
  ))
}
