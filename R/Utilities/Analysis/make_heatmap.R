#' Create heatmap with optional feature selection (ANOVA / variance / MAD / t-test)
#' and Stage annotation. Returns plot object for patchwork.
#'
#' @param data                 Data frame with ID, stage_bin, then feature columns
#' @param variant_colors       Named color vector for Variant (names = levels)
#' @param top_features         NULL (default: show all features). If numeric >0, keep top N by `feature_selector`.
#' @param feature_selector     One of c("none","anova","variance","mad","ttest"). Default "none".
#' @param variant_levels       Factor order for Variant (default c("PTC","FV-PTC","FTC"))
#' @param stage_colors         Named color vector for stage_bin; include "Early" and "Advanced"
#' @return List with plot object, M, Mz, hc_cols, ann_col, ann_colors, etc.
#' @export
make_heatmap <- function(
    data,
    variant_colors = c("PTC" = "#DF8D0A", "FV-PTC" = "#23744E", "FTC" = "#194992"),
    top_features = NULL,
    feature_selector = c("none", "anova", "variance", "mad", "ttest"),
    variant_levels = c("PTC", "FV-PTC", "FTC"),
    stage_colors = c("Early" = "#113d6a", "Advanced" = "#800017")) {
  feature_selector <- match.arg(feature_selector)

  # ---- Checks ----
  stopifnot(all(c("ID", "stage_bin") %in% names(data)))

  # Derive Variant from ID
  dat <- data |>
    dplyr::mutate(
      Variant = dplyr::case_when(
        grepl("^FVPTC", ID) ~ "FV-PTC",
        grepl("^F[0-9]", ID) ~ "FTC",
        grepl("^P[0-9]", ID) ~ "PTC",
        TRUE ~ NA_character_
      )
    )
  
  # Coerce Variant to factor in desired order
  dat$Variant <- factor(dat$Variant, levels = variant_levels)

  # Build matrix: samples x features (drop ID/stage_bin/Variant)
  drop_cols <- c("ID", "stage_bin", "Variant")
  X <- as.matrix(dplyr::select(dat, -dplyr::all_of(drop_cols)))
  stopifnot(all(vapply(as.data.frame(X), is.numeric, TRUE)))

  # Sample IDs (use make.names for uniqueness + valid rownames)
  sample_ids <- make.names(dat$ID, unique = TRUE)
  rownames(X) <- sample_ids

  # Group factor for ANOVA ranking (Variant)
  group_variant <- dat$Variant
  
  # Stage group for t-test ranking
  group_stage <- dat$stage_bin

  # ---- Optional feature ranking/selection ----
  if (!is.null(top_features) && is.numeric(top_features) && top_features > 0 &&
    feature_selector != "none") {
    top_n <- min(top_features, ncol(X))

    rank_idx <- switch(feature_selector,
      "anova" = {
        pvals <- apply(X, 2, function(x) {
          fit <- aov(x ~ group_variant)
          summary(fit)[[1]][["Pr(>F)"]][1]
        })
        order(pvals, na.last = TRUE) # ascending p
      },
      "ttest" = {
        pvals <- apply(X, 2, function(x) {
          tryCatch({
            t.test(x ~ group_stage)$p.value
          }, error = function(e) NA_real_)
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

  # ---- Column annotation (aligned to columns of M) ----
  # ---- Column annotation (aligned to columns of M) ----
  # Build annotation in REVERSE order since pheatmap displays first column at bottom
  # Start with variant annotation (will be on bottom - first column)
  ann_col <- data.frame(Variant = dat$Variant, row.names = sample_ids)

  # Add Stage annotation (will be on top - second column)
  # Convert stage_bin to factor with canonical order
  stage_factor <- factor(as.character(dat$stage_bin), levels = c("Early", "Advanced"))
  ann_col$Stage <- stage_factor

  # Reorder rows of ann_col to match M's columns
  ann_col <- ann_col[colnames(M), , drop = FALSE]

  # ---- Annotation color lists ----
  # Build colors in same order as ann_col data frame (pheatmap displays in reverse)
  ann_colors <- list()
  
  # Add variant colors first (matches first column - will display at bottom)
  ann_colors$Variant <- variant_colors

  # Add stage colors (matches second column - will display at top)
  ann_colors$Stage <- stage_colors

  # ---- Heatmap (for screen) ----
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
    na_col = "#DDDDDD",
    legend_labels = "Z-Score"
  )

  # Create heatmap plot object for patchwork
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
    silent = TRUE,  # Prevents auto-display
    legend_labels = "Z-Score"
  )

  list(
    M = M,
    Mz = Mz,
    hc_cols = hc_cols,
    sample_ids = sample_ids,
    group_variant = group_variant,
    group_stage = group_stage,
    feature_selector = feature_selector,
    top_features = top_features,
    ann_col = ann_col,
    ann_colors = ann_colors,
    heatmap_plot = heatmap_plot  # Plot object for patchwork
  )
}
