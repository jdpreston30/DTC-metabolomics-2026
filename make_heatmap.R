#' Create heatmap from feature table with annotation
#'
#' @param data                 Feature table (tibble or data frame) with ID column, annotation column (e.g., stage_bin), then feature columns
#' @param id_col              Name of the ID column (default "ID")
#' @param annotation_col      Name of the annotation column (default "stage_bin")
#' @param annotation_colors   Named color vector for annotation levels
#' @param top_features        NULL (default: show all features). If numeric >0, keep top N by `feature_selector`.
#' @param feature_selector    One of c("none","ttest","anova","variance","mad"). Default "none".
#'
#' @return List with plot object, M, Mz, hc_cols, ann_col, ann_colors, etc.
#' @export
make_heatmap <- function(
    data,
    id_col = "ID",
    annotation_col = "stage_bin",
    annotation_colors = NULL,
    top_features = NULL,
    feature_selector = c("none", "ttest", "anova", "variance", "mad")) {
  feature_selector <- match.arg(feature_selector)

  # ---- Checks ----
  stopifnot(id_col %in% names(data))
  stopifnot(annotation_col %in% names(data))

  # Convert to data frame for easier manipulation
  dat <- as.data.frame(data)

  # Build matrix: samples x features (drop ID and annotation columns)
  drop_cols <- c(id_col, annotation_col)
  X <- as.matrix(dplyr::select(dat, -dplyr::all_of(drop_cols)))
  stopifnot(all(vapply(as.data.frame(X), is.numeric, TRUE)))

  # Sample IDs (use make.names for uniqueness + valid rownames)
  sample_ids <- make.names(dat[[id_col]], unique = TRUE)
  rownames(X) <- sample_ids

  # Group factor for ANOVA ranking
  group <- factor(dat[[annotation_col]])

  # ---- Optional feature ranking/selection ----
  if (!is.null(top_features) && is.numeric(top_features) && top_features > 0 &&
    feature_selector != "none") {
    top_n <- min(top_features, ncol(X))

    rank_idx <- switch(feature_selector,
      "ttest" = {
        # T-test for two-group comparison
        group_levels <- levels(group)
        if (length(group_levels) != 2) {
          stop("T-test feature selector requires exactly 2 groups")
        }
        pvals <- apply(X, 2, function(x) {
          g1 <- x[group == group_levels[1]]
          g2 <- x[group == group_levels[2]]
          tryCatch({
            t.test(g1, g2)$p.value
          }, error = function(e) NA_real_)
        })
        order(pvals, na.last = TRUE) # ascending p
      },
      "anova" = {
        pvals <- apply(X, 2, function(x) {
          fit <- aov(x ~ group)
          summary(fit)[[1]][["Pr(>F)"]][1]
        })
        order(pvals, na.last = TRUE) # ascending p
      },
      "variance" = {
        v <- apply(X, 2, stats::var, na.rm = TRUE)
        order(v, decreasing = TRUE, na.last = NA)
      },
      "mad" = {
        m <- apply(X, 2, stats::mad, na.rm = TRUE)
        order(m, decreasing = TRUE, na.last = NA)
      }
    )

    X <- X[, head(rank_idx, top_n), drop = FALSE]
  }

  # Drop zero-variance columns (after selection)
  nzv <- apply(X, 2, sd, na.rm = TRUE) > 0
  if (!all(nzv)) X <- X[, nzv, drop = FALSE]
  if (!ncol(X)) stop("No features remain after selection/variance filtering.")

  # Heatmap matrix: features x samples
  M <- t(X)

  # ---- Sample clustering ----
  Mz <- t(scale(t(M), center = TRUE, scale = TRUE))
  Mz[is.na(Mz)] <- 0
  d_cols <- dist(t(Mz), method = "euclidean")
  hc_cols <- hclust(d_cols, method = "complete")

  # ---- Column annotation ----
  ann_col <- data.frame(
    Annotation = factor(dat[[annotation_col]]),
    row.names = sample_ids
  )
  colnames(ann_col) <- annotation_col
  
  # Reorder to match M's columns
  ann_col <- ann_col[colnames(M), , drop = FALSE]

  # ---- Annotation colors ----
  ann_colors <- list()
  
  if (is.null(annotation_colors)) {
    # Generate default colors based on number of levels
    levels_n <- nlevels(ann_col[[annotation_col]])
    default_colors <- RColorBrewer::brewer.pal(
      max(3, min(levels_n, 9)), 
      "Set2"
    )[1:levels_n]
    annotation_colors <- setNames(default_colors, levels(ann_col[[annotation_col]]))
  }
  
  ann_colors[[annotation_col]] <- annotation_colors

  # ---- Heatmap ----
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
    fontsize = 8,
    na_col = "#DDDDDD",
    silent = TRUE,
    legend_labels = "Z-Score"
  )

  list(
    M = M,
    Mz = Mz,
    hc_cols = hc_cols,
    sample_ids = sample_ids,
    group = group,
    feature_selector = feature_selector,
    top_features = top_features,
    ann_col = ann_col,
    ann_colors = ann_colors,
    heatmap_plot = heatmap_plot
  )
}
