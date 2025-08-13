#' Create heatmap with optional feature selection (ANOVA / variance / MAD)
#' and optional T_stage annotation. Adds optional PNG export.
#'
#' @param data             Data frame with Patient_ID, Variant, optional T_stage, then feature columns
#' @param variant_colors   Named color vector for Variant (names = levels)
#' @param top_features     NULL (default: show all features). If numeric >0, keep top N by `feature_selector`.
#' @param feature_selector One of c("none","anova","variance","mad"). Default "none".
#' @param variant_levels   Factor order for Variant (default c("PTC","FV-PTC","FTC"))
#' @param n_clades         Number of clades to extract from sample HCA (default 2)
#' @param output_filename  Optional filename (saved to "Outputs/") for SVG
#' @param annotate_t_stage Logical; if TRUE and T_stage exists, add as column annotation
#' @param T_stage_colors   Optional named color vector for T_stage; include "T1","T2","T3","T4". "Unknown" added if missing.
#' @param output_png_filename Optional PNG filename (saved to "Outputs/"). If NULL, no PNG is written.
#' @param png_width        PNG width in pixels (default 2000)
#' @param png_height       PNG height in pixels (default 2000)
#' @param png_res          PNG resolution (dpi-ish; default 220)
#' @param png_bg           PNG background color (default "white")
#' @return List with M, Mz, hc_cols, clade_df, clade_lists, ann_col, ann_colors, etc.
#' @export
make_heatmap <- function(
    data,
    variant_colors,
    top_features = NULL,
    feature_selector = c("none", "anova", "variance", "mad"),
    variant_levels = c("PTC", "FV-PTC", "FTC"),
    n_clades = 2,
    output_filename = NULL, # SVG
    annotate_t_stage = FALSE,
    T_stage_colors = NULL,
    # --- NEW: PNG options ---
    output_png_filename = NULL, # e.g., "heatmap.png"
    png_width = 2000,
    png_height = 2000,
    png_res = 220,
    png_bg = "white") {
  feature_selector <- match.arg(feature_selector)

  # ---- Checks ----
  stopifnot(all(c("Patient_ID", "Variant") %in% names(data)))
  has_t <- "T_stage" %in% names(data)

  # Keep ID/Variant/(optional)T_stage up front
  dat <- dplyr::select(
    data,
    dplyr::any_of(c("Patient_ID", "Variant", if (has_t) "T_stage")),
    dplyr::everything()
  )

  # Coerce Variant to factor in desired order
  dat$Variant <- factor(dat$Variant, levels = variant_levels)

  # Build matrix: samples x features (drop ID/Variant/T_stage)
  drop_cols <- c("Patient_ID", "Variant", if (has_t) "T_stage")
  X <- as.matrix(dplyr::select(dat, -dplyr::all_of(drop_cols)))
  stopifnot(all(vapply(as.data.frame(X), is.numeric, TRUE)))

  # Sample IDs (use make.names for uniqueness + valid rownames)
  sample_ids <- make.names(dat$Patient_ID, unique = TRUE)
  rownames(X) <- sample_ids

  # Group factor for ANOVA ranking
  group <- dat$Variant

  # ---- Optional feature ranking/selection ----
  if (!is.null(top_features) && is.numeric(top_features) && top_features > 0 &&
    feature_selector != "none") {
    top_n <- min(top_features, ncol(X))

    rank_idx <- switch(feature_selector,
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

  # ---- Column annotation (aligned to columns of M) ----
  ann_col <- data.frame(Variant = dat$Variant, row.names = sample_ids)

  if (annotate_t_stage && has_t) {
    # Turn T_stage into character, map NA -> "Unknown"
    t_raw <- as.character(dat$T_stage)
    t_raw[is.na(t_raw)] <- "Unknown"

    # If values already start with "T" (e.g., "T3"), keep them;
    # if they look numeric (e.g., "3"), prepend "T".
    needs_T <- grepl("^[0-9]+$", t_raw)
    t_clean <- ifelse(needs_T, paste0("T", t_raw), t_raw)

    # Factor with canonical order (include Unknown)
    t_lvls <- c("T1", "T2", "T3", "T4", "Unknown")
    ann_col$T_stage <- factor(t_clean, levels = t_lvls)
  }

  # Reorder rows of ann_col to match M's columns
  ann_col <- ann_col[colnames(M), , drop = FALSE]

  # ---- Annotation color lists ----
  ann_colors <- list(Variant = variant_colors)

  if (annotate_t_stage && has_t) {
    # Seed T-stage colors if not provided
    if (is.null(T_stage_colors)) T_stage_colors <- c()
    # Ensure Unknown exists
    if (!("Unknown" %in% names(T_stage_colors))) {
      T_stage_colors <- c(T_stage_colors, Unknown = "#BDBDBD")
    }
    # Restrict to levels that appear; fill any missing with greys
    needed <- levels(ann_col$T_stage)
    missing <- setdiff(needed, names(T_stage_colors))
    if (length(missing)) {
      fill_cols <- grDevices::gray.colors(length(missing), start = 0.3, end = 0.7)
      names(fill_cols) <- missing
      T_stage_colors <- c(T_stage_colors, fill_cols)
    }
    ann_colors$T_stage <- T_stage_colors[needed]
  }

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
    na_col = "#DDDDDD"
  )

  # ---- Optional save to SVG ----
  if (!is.null(output_filename)) {
    dir.create("Outputs", showWarnings = FALSE)
    grDevices::svg(file.path("Outputs", output_filename), width = 8, height = 8, pointsize = 12)
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
    grDevices::dev.off()
    cat("Heatmap saved to:", file.path("Outputs", output_filename), "\n")
  }

  # ---- NEW: Optional save to PNG ----
  if (!is.null(output_png_filename)) {
    dir.create("Outputs", showWarnings = FALSE)
    grDevices::png(
      filename = file.path("Outputs", output_png_filename),
      width = png_width, height = png_height,
      res = png_res, units = "px", bg = png_bg
    )
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
    grDevices::dev.off()
    cat("Heatmap PNG saved to:", file.path("Outputs", output_png_filename), "\n")
  }

  # ---- Sample clustering & clades ----
  Mz <- t(scale(t(M), center = TRUE, scale = TRUE))
  Mz[is.na(Mz)] <- 0
  d_cols <- dist(t(Mz), method = "euclidean")
  hc_cols <- hclust(d_cols, method = "complete")

  # Map back to original Patient_ID
  id_map <- setNames(dat$Patient_ID, sample_ids)
  clades <- stats::cutree(hc_cols, k = n_clades)
  ids_ordered_clean <- names(clades)[hc_cols$order]
  ids_ordered_orig <- unname(id_map[ids_ordered_clean])

  clade_df <- tibble::tibble(
    Patient_ID = ids_ordered_orig,
    Clade      = unname(clades[ids_ordered_clean])
  )

  clade_lists <- lapply(seq_len(n_clades), function(i) {
    clade_df %>%
      dplyr::filter(Clade == i) %>%
      dplyr::pull(Patient_ID)
  })
  names(clade_lists) <- paste0("clade", seq_len(n_clades), "_ids")

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
    clade_df = clade_df,
    clades = clades,
    clade_lists = clade_lists
  )
}
