#' Create heatmap with optional feature selection (ANOVA / variance / MAD / t-test)
#' and Stage annotation. Returns plot object for patchwork.
#'
#' @param data                 Data frame with ID, stage_bin, then feature columns
#' @param variant_colors       Named color vector for Variant (names = levels)
#' @param top_features         NULL (default: show all features). If numeric >0, keep top N by `feature_selector`.
#' @param feature_selector     One of c("none","anova","variance","mad","ttest"). Default "none".
#' @param variant_levels       Factor order for Variant (default c("PTC","FV-PTC","FTC"))
#' @param stage_colors         Named color vector for stage_bin; include "Early" and "Advanced"
#' @param print_preview        Logical. If TRUE, saves a PNG preview of the heatmap. Default FALSE.
#' @param print_scale           Logical. If TRUE, shows the color scale legend. Default TRUE.
#' @param alt_variables        NULL (default: use Variant+Stage). If "T_stage_comp", use T_stage_comp annotation instead.
#' @param annotation_height    Numeric. Total height of annotation bars. Default NA (auto-size).
#' @return List with plot object, M, Mz, hc_cols, ann_col, ann_colors, etc.
#' @export
make_heatmap <- function(
    data,
    variant_colors = c("PTC" = "#DF8D0A", "FV-PTC" = "#23744E", "FTC" = "#194992"),
    top_features = NULL,
    feature_selector = c("none", "anova", "variance", "mad", "ttest"),
    variant_levels = c("PTC", "FV-PTC", "FTC"),
    stage_colors = c("Early" = "#113d6a", "Advanced" = "#800017"),
    print_preview = FALSE,
    print_scale = TRUE,
    alt_variables = NULL,
    annotation_height = NA) {
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
  
  # If using alt_variables, also drop that column from features
  if (!is.null(alt_variables) && alt_variables == "T_stage_comp") {
    drop_cols <- c(drop_cols, "T_stage_comp")
  }
  
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
  # Build annotation in REVERSE order since pheatmap displays first column at bottom
  # Start with variant annotation (will be on bottom - first column)
  ann_col <- data.frame(Variant = dat$Variant, row.names = sample_ids)

  # ---- Annotation color lists ----
  # Build colors in same order as ann_col data frame (pheatmap displays in reverse)
  ann_colors <- list()
  
  # Add variant colors first (matches first column - will display at bottom)
  ann_colors$Variant <- variant_colors

  # If using alt_variables, add T_stage_comp annotation next (will be middle)
  if (!is.null(alt_variables) && alt_variables == "T_stage_comp") {
    # Verify T_stage_comp column exists
    if (!"T_stage_comp" %in% names(dat)) {
      stop("T_stage_comp column not found in data. Required when alt_variables = 'T_stage_comp'")
    }
    
    # Add T_stage_comp annotation with custom label (second column - will display in middle)
    ann_col$`T Category` <- factor(dat$T_stage_comp, levels = c("T1", "T2", "T3", "T4"), ordered = TRUE)
    
    # Add sequential heat color palette (cool to hot: T1 -> T4)
    # Using yellow-orange-red progression distinct from Stage colors
    ann_colors$`T Category` <- c(
      "T1" = "#ffffcc",  # Light yellow (coolest)
      "T2" = "#ffeda0",  # Pale yellow
      "T3" = "#feb24c",  # Orange
      "T4" = "#f03b20"   # Bright red (hottest)
    )
  }

  # Add Stage annotation last (will be on top)
  # Convert stage_bin to factor with canonical order
  stage_factor <- factor(as.character(dat$stage_bin), levels = c("Early", "Advanced"))
  ann_col$Stage <- stage_factor
  
  # Add stage colors
  ann_colors$Stage <- stage_colors

  # Reorder rows of ann_col to match M's columns
  ann_col <- ann_col[colnames(M), , drop = FALSE]

  # ---- Create heatmap plot object (single call with silent=TRUE) ----
  # Define color palette
  heatmap_colors <- colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(255)
  
  # Build arguments list
  pheatmap_args <- list(
    mat = M,
    scale = "row",
    color = heatmap_colors,
    clustering_distance_rows = "euclidean",
    clustering_distance_cols = "euclidean",
    clustering_method = "complete",
    annotation_col = ann_col,
    annotation_colors = ann_colors,
    annotation_names_col = FALSE,
    annotation_legend = TRUE,
    show_rownames = FALSE,
    show_colnames = FALSE,
    fontsize = 6.55,
    na_col = "#DDDDDD",
    silent = TRUE,
    legend = print_scale,
    legend_labels = "Z-Score"
  )
  
  if (!is.na(annotation_height)) {
    pheatmap_args$annotation_height <- annotation_height
  }
  
  heatmap_plot <- do.call(pheatmap::pheatmap, pheatmap_args)
  
  # Hide all legend titles recursively
  hide_titles <- function(grob) {
    if (inherits(grob, "titleGrob")) {
      # This is a title grob - make it transparent
      if (!is.null(grob$gp)) {
        grob$gp$col <- "transparent"
      } else {
        grob$gp <- grid::gpar(col = "transparent")
      }
    }
    if (inherits(grob, "text")) {
      # Check if this is in a legend title position
      if (!is.null(grob$label) && length(grob$label) > 0) {
        label_val <- as.character(grob$label)
        if (length(label_val) == 1 && label_val %in% c("Stage", "T Category", "Variant")) {
          if (!is.null(grob$gp)) {
            grob$gp$col <- "transparent"
          } else {
            grob$gp <- grid::gpar(col = "transparent")
          }
        }
      }
    }
    if (inherits(grob, "gTree") && !is.null(grob$children)) {
      grob$children <- lapply(grob$children, hide_titles)
    }
    return(grob)
  }
  
  heatmap_plot$gtable$grobs <- lapply(heatmap_plot$gtable$grobs, hide_titles)
  
  # Store legend parameters for custom legend creation
  legend_params <- list(
    colors = heatmap_colors,
    limits = c(-3, 3),  # Typical z-score range
    breaks = seq(-3, 3, by = 1),
    labels = c("-3", "-2", "-1", "0", "1", "2", "3"),
    title = "Z-Score",
    fontsize = 8  # Match heatmap fontsize
  )

  # Optional print preview
  if (print_preview) {
    preview_filename <- paste0(feature_selector, "_", 
                              ifelse(!is.null(top_features), paste0("top", top_features), "all"), 
                              ".png")
    print_to_png(
      plot = heatmap_plot,
      filename = preview_filename,
      dpi = 300,
      width = 6,
      height = 6
    )
  }

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
    heatmap_plot = heatmap_plot,  # Full pheatmap object
    legend_params = legend_params  # Parameters for custom legend
  )
}
