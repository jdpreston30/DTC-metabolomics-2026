# ===== File: Analysis/build_enrichment_network.R =====
#' Build an Enrichment Network Plot
#'
#' This function creates an EnrichNet-style network plot from enrichment results,
#' where nodes represent pathways and edges represent similarity based on shared compounds.
#' Node colors correspond to p-values, and node sizes correspond to enrichment factors.
#' The plot includes customizable legends for enrichment factors and p-values,
#' and supports saving the output to a file.
#'
#' @param enrich_df A data frame containing enrichment results. Must include `pathway_ID`, `pathway_name`, `enrichment_factor`, and either `p_value` or `neg_log_p`.
#' @param edge_thresh Numeric threshold for edge inclusion based on Jaccard similarity of compounds (default 0.10).
#' @param prefer_hsa Logical indicating whether to prefer human (hsa) KEGG pathway IDs when fetching compounds (default TRUE).
#' @param term2compound_override Optional data frame to override pathway-to-compound mappings. Must have columns `pathway_ID` and `compound_id` (default NULL).
#' @param save_path Optional file path to save the plot (default NULL).
#' @param plot_title Optional title for the plot (default NULL).
#' @param width Numeric width of saved plot in units (default 8).
#' @param height Numeric height of saved plot in units (default 6).
#' @param dpi Numeric resolution for saved plot (default 300).
#' @param units Units for width and height when saving (default "in").
#' @param bg Background color for saved plot (default "white").
#' @param seed Integer seed for reproducible layout (default 123).
#' @param layout Layout algorithm for graph plotting (default "fr").
#' @param p_limits Numeric vector of length 2 specifying limits for the p-value color scale (default c(0.01, 0.05)).
#' @param show_enrichment Logical to toggle the enrichment factor legend (default TRUE).
#' @param show_pvalue Logical to toggle the p-value legend (default TRUE).
#'
#' @return A list containing:
#' \describe{
#'   \item{plot}{The ggplot2 object of the enrichment network.}
#'   \item{graph}{The igraph object representing the network.}
#'   \item{nodes}{Data frame of node attributes.}
#'   \item{edges}{Data frame of edge attributes.}
#'   \item{term2compound}{Data frame mapping pathways to compounds.}
#' }
#'
#' @examples
#' \dontrun{
#'   # Assuming enrich_df is a data frame with required columns
#'   result <- build_enrichment_network(enrich_df)
#'   print(result$plot)
#' }
build_enrichment_network <- function(
    enrich_df,
    edge_thresh = 0.10,
    prefer_hsa = TRUE,
    term2compound_override = NULL,
    save_path = NULL,
    plot_title = NULL,
    width = 8, height = 6, dpi = 300, units = "in", bg = "white",
    seed = 123, layout = "fr",
    p_limits = c(0.01, 0.05),
    show_enrichment = TRUE,   # NEW: toggle enrichment factor legend
    show_pvalue     = TRUE    # NEW: toggle p-value legend
) {
  pkgs <- c("dplyr", "tidyr", "purrr", "stringr", "igraph", "ggraph", "scales", "ggplot2", "KEGGREST", "memoise", "grid")
  missing <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
  if (length(missing)) stop("Install missing packages: ", paste(missing, collapse = ", "))

  # ---- prep input: ensure we have p_value (derive from neg_log_p if needed) ----
  if (!"p_value" %in% names(enrich_df)) {
    if ("neg_log_p" %in% names(enrich_df)) {
      enrich_df <- dplyr::mutate(enrich_df, p_value = 10^(-.data$neg_log_p))
    } else {
      stop("enrich_df must contain either 'p_value' or 'neg_log_p'.")
    }
  }

  # ---- KEGG membership (same as before) ----
  kegg_link_m <- memoise::memoise(function(kind, x) KEGGREST::keggLink(kind, x))
  to_hsa <- function(pid) stringr::str_replace(pid, "^map", "hsa")
  get_pathway_compounds <- function(pid) {
    candidates <- if (prefer_hsa) unique(c(to_hsa(pid), pid)) else pid
    for (cand in candidates) {
      res <- tryCatch(kegg_link_m("cpd", paste0("path:", cand)), error = function(e) NULL)
      if (!is.null(res) && length(res) > 0) {
        return(list(path = cand, compounds = sub("^cpd:", "", unname(res))))
      }
    }
    list(path = if (prefer_hsa) to_hsa(pid) else pid, compounds = character(0))
  }

  if (is.null(term2compound_override)) {
    t2c <- enrich_df |>
      dplyr::distinct(pathway_ID) |>
      dplyr::mutate(fetch = purrr::map(pathway_ID, get_pathway_compounds)) |>
      dplyr::mutate(
        final_id = purrr::map_chr(fetch, "path"),
        compounds = purrr::map(fetch, "compounds")
      ) |>
      dplyr::select(pathway_ID, final_id, compounds) |>
      tidyr::unnest_longer(compounds, values_to = "compound_id") |>
      dplyr::mutate(compound_id = as.character(compound_id)) |>
      dplyr::filter(!is.na(compound_id), nzchar(compound_id)) |>
      dplyr::mutate(pathway_ID = final_id) |>
      dplyr::select(-final_id) |>
      dplyr::distinct()
  } else {
    t2c <- term2compound_override |>
      dplyr::select(pathway_ID, compound_id) |>
      dplyr::distinct()
  }

  t2c <- t2c |> dplyr::semi_join(enrich_df, by = "pathway_ID")

  sets <- t2c |>
    dplyr::group_by(pathway_ID) |>
    dplyr::summarise(members = list(unique(compound_id)), .groups = "drop")

  pairs <- tidyr::crossing(a = sets$pathway_ID, b = sets$pathway_ID) |>
    dplyr::filter(a < b) |>
    dplyr::left_join(sets |> dplyr::rename(members_a = members), by = c("a" = "pathway_ID")) |>
    dplyr::left_join(sets |> dplyr::rename(members_b = members), by = c("b" = "pathway_ID")) |>
    dplyr::mutate(
      inter   = purrr::map2_int(members_a, members_b, ~ length(intersect(.x, .y))),
      union_n = purrr::map2_int(members_a, members_b, ~ length(union(.x, .y))),
      jaccard = dplyr::if_else(union_n > 0, inter / union_n, 0)
    )

  edges <- pairs |>
    dplyr::filter(jaccard >= edge_thresh) |>
    dplyr::transmute(from = a, to = b, weight = jaccard)

  nodes <- enrich_df |>
    dplyr::rename(
      kegg_id = pathway_ID, Name = pathway_name,
      size_val = enrichment_factor, p_val = p_value
    ) |>
    dplyr::distinct() |>
    dplyr::mutate(
      Name = stringr::str_replace_all(
        Name,
        stringr::regex("\\b(Degradation|Metabolism|Biosynthesis)\\b", ignore_case = TRUE),
        "\n\\1"
      ),
      Name = stringr::str_replace_all(Name, "\\s*\\n\\s*", "\n") # tidy spaces around break
    )

  g <- igraph::graph_from_data_frame(
    d = edges,
    vertices = nodes |> dplyr::select(kegg_id, Name, size_val, p_val),
    directed = FALSE
  )

  # ---- Compute layout coordinates ----
  set.seed(seed)
  coords <- switch(
    layout,
    fr = igraph::layout_with_fr(g),
    kk = igraph::layout_with_kk(g),
    lgl = igraph::layout_with_lgl(g),
    circle = igraph::layout_in_circle(g),
    grid = igraph::layout_on_grid(g),
    stop("Unsupported layout: ", layout)
  )
  coords_df <- as.data.frame(coords)
  colnames(coords_df) <- c("x", "y")
  coords_df$name <- igraph::V(g)$name

  # ---- PLOT ----
  p <- ggraph::ggraph(g, layout = coords_df) +
    ggraph::geom_edge_link(aes(width = weight), alpha = 0.5, colour = "grey50") +

    # Nodes
    ggraph::geom_node_point(aes(size = size_val, fill = p_val),
      shape = 21, stroke = 0.6, colour = "black"
    ) +

    # Labels
    ggraph::geom_node_text(aes(label = Name),
      size = 4, colour = "black",
      family = "Arial", fontface = "bold",
      vjust = 2.2
    ) +

    ggraph::scale_edge_width(range = c(0.2, 2), guide = "none") +
    theme_void(base_family = "Arial") +
    theme(
      legend.position = "right",
      legend.margin   = margin(l = 40, r = 0, t = 0, b = 0),
      legend.text     = element_text(size = 13, face = "bold", family = "Arial"),
      legend.title    = element_text(size = 14, face = "bold", family = "Arial"),
      text            = element_text(color = "black", family = "Arial")
    ) +
    ggplot2::coord_equal(clip = "off") +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = 0.07)) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = 0.07)) +
    {if (!is.null(plot_title)) ggplot2::labs(title = plot_title)} +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        size   = 25,
        face   = "bold",
        family = "Arial",
        hjust  = 0.5,
        vjust  = 5
      ),
      plot.margin = ggplot2::margin(5, 5, 5, 5)
    )

  # --- Conditional legends ---
  p <- p +
    scale_size_continuous(
      range = c(3, 30),
      limits = c(0, 5),
      breaks = c(4, 3, 2, 1),
      name = "\nEnrichment factor",
      guide = if (show_enrichment) {
        guide_legend(
          reverse = TRUE,
          override.aes = list(fill = "black", colour = "black")
        )
      } else {
        "none"
      }
    )

  p <- p +
    scale_fill_gradient(
      low = "#0a2256", high = "#c3dbe9",
      limits = p_limits, oob = scales::squish,
      name = "p-value\n",
      guide = if (show_pvalue) {
        guide_colorbar(
          reverse   = TRUE,
          barheight = grid::unit(8, "cm"),
          barwidth  = grid::unit(1.3, "cm")
        )
      } else {
        "none"
      }
    )

  # --- Saving ---
  if (!is.null(save_path)) {
    ggplot2::ggsave(
      filename = save_path, plot = p,
      width = width, height = height,
      dpi = dpi, units = units, bg = bg
    )
  }

  list(plot = p, graph = g, nodes = nodes, edges = edges, term2compound = t2c)
}

# ===== File: Analysis/get_permanova.R =====
#' Run individual variable PERMANOVA analysis
#'
#' This function performs PERMANOVA (Permutational Multivariate Analysis of Variance)
#' for a single variable against a distance matrix. It's designed to be used with
#' lapply to test multiple variables individually.
#'
#' @param varname Character string specifying the variable name to test
#' @param feature_data Data frame containing the feature/metabolite data (samples x features)
#' @param meta_data Data frame containing the metadata with all variables
#' @param ctrl Permutation control object from permute::how()
#' @param distance_method Character string specifying distance method for vegan::vegdist. Default "euclidean"
#' @param seed Integer for random seed to ensure reproducibility. Default 2025
#'
#' @return A tibble with columns:
#'   \item{Variable}{Character, the variable name tested}
#'   \item{R2}{Numeric, the R-squared value from PERMANOVA}
#'   \item{p_value}{Numeric, the p-value from permutation test}
#'
#' @details
#' The function:
#' \itemize{
#'   \item Sets a seed for reproducibility
#'   \item Removes samples with missing values for the test variable
#'   \item Skips variables with insufficient factor levels (< 2)
#'   \item Calculates distance matrix using specified method
#'   \item Runs adonis2 PERMANOVA test
#'   \item Returns results in tidy format
#' }
#'
#' @examples
#' \dontrun{
#' # Set up data
#' features <- UFT_data[, feature_cols]
#' metadata <- UFT_data[, c("Sex", "Age", "Variant")]
#' ctrl <- permute::how(nperm = 999)
#' 
#' # Test single variable
#' result <- get_permanova("Sex", features, metadata, ctrl)
#' 
#' # Test multiple variables
#' variables <- c("Sex", "Age", "Variant")
#' results <- bind_rows(lapply(variables, get_permanova, 
#'                            feature_data = features, 
#'                            meta_data = metadata, 
#'                            ctrl = ctrl))
#' }
#'
#' @export
get_permanova <- function(varname, 
                         feature_data, 
                         meta_data, 
                         ctrl, 
                         distance_method = "euclidean", 
                         seed = 2025) {
  
  # Set seed for reproducibility
  set.seed(seed)
  
  # Check that variable exists in metadata
  if (!varname %in% names(meta_data)) {
    stop(paste("Variable", varname, "not found in meta_data"))
  }
  
  # Remove samples with missing values for the test variable
  keep <- complete.cases(meta_data[[varname]])
  sub_meta <- droplevels(meta_data[keep, , drop = FALSE])
  sub_features <- feature_data[keep, , drop = FALSE]
  
  # Skip variables with insufficient levels for factors
  if (is.factor(sub_meta[[varname]]) &&
      nlevels(droplevels(sub_meta[[varname]])) < 2) {
    return(tibble::tibble(
      Variable = varname, 
      R2 = NA_real_, 
      p_value = NA_real_
    ))
  }
  
  # Calculate distance matrix
  d <- vegan::vegdist(sub_features, method = distance_method)
  
  # Run PERMANOVA
  formula_str <- paste("d ~ `", varname, "`", sep = "")
  res <- vegan::adonis2(
    as.formula(formula_str), 
    data = sub_meta, 
    permutations = ctrl
  )
  
  # Return results in tidy format
  tibble::tibble(
    Variable = varname,
    R2 = as.numeric(res$R2[1]),
    p_value = as.numeric(res$`Pr(>F)`[1])
  )
}

# ===== File: Analysis/make_heatmap.R =====
#' Create heatmap with optional feature selection (ANOVA / variance / MAD)
#' and optional T_stage annotation. Returns plot object for patchwork.
#'
#' @param data                 Data frame with Patient_ID, Variant, optional T_stage, then feature columns
#' @param variant_colors       Named color vector for Variant (names = levels)
#' @param top_features         NULL (default: show all features). If numeric >0, keep top N by `feature_selector`.
#' @param feature_selector     One of c("none","anova","variance","mad"). Default "none".
#' @param variant_levels       Factor order for Variant (default c("PTC","FV-PTC","FTC"))
#' @param n_clades             Number of clades to extract from sample HCA (default 2)
#' @param annotate_t_stage     Logical; if TRUE and T_stage exists, add as column annotation
#' @param T_stage_colors       Optional named color vector for T_stage; include "T1","T2","T3","T4". "Unknown" added if missing.
#' @param cluster_colors       Named color vector for cluster annotation (e.g., c("Cluster 1" = "#color1", "Cluster 2" = "#color2"))
#' @return List with plot object, M, Mz, hc_cols, clade_df, clade_lists, ann_col, ann_colors, etc.
#' @export
make_heatmap <- function(
    data,
    variant_colors = c("PTC" = "#DF8D0A", "FV-PTC" = "#23744E", "FTC" = "#194992"),
    top_features = NULL,
    feature_selector = c("none", "anova", "variance", "mad"),
    variant_levels = c("PTC", "FV-PTC", "FTC"),
    n_clades = 2,
    annotate_t_stage = FALSE,
    T_stage_colors = c("T1-T2" = "#dfba37", "T3-T4" = "#72061c"),
    cluster_colors = c("Cluster 1" = "#94001E", "Cluster 2" = "#03507D")) {
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

  # ---- Sample clustering & clades ----
  Mz <- t(scale(t(M), center = TRUE, scale = TRUE))
  Mz[is.na(Mz)] <- 0
  d_cols <- dist(t(Mz), method = "euclidean")
  hc_cols <- hclust(d_cols, method = "complete")

  # Map back to original Patient_ID
  id_map <- setNames(dat$Patient_ID, sample_ids)
  clades_raw <- stats::cutree(hc_cols, k = n_clades)
  
  # Assign clusters based on dendrogram order: leftmost = Cluster 1, rightmost = Cluster 2
  # Get the order of samples from the dendrogram
  ordered_samples <- names(clades_raw)[hc_cols$order]
  
  # Find which raw cluster appears first (leftmost) in the dendrogram order
  first_cluster_raw <- clades_raw[ordered_samples[1]]
  
  # Assign final cluster numbers: leftmost cluster becomes Cluster 1, rightmost becomes Cluster 2
  clades <- ifelse(clades_raw == first_cluster_raw, 1, 2)
  
  ids_ordered_clean <- names(clades)[hc_cols$order]
  ids_ordered_orig <- unname(id_map[ids_ordered_clean])

  cluster_df <- tibble::tibble(
    Patient_ID = ids_ordered_orig,
    Cluster    = unname(clades[ids_ordered_clean])
  )

  cluster_lists <- lapply(seq_len(n_clades), function(i) {
    cluster_df %>%
      dplyr::filter(Cluster == i) %>%
      dplyr::pull(Patient_ID)
  })
  names(cluster_lists) <- paste0("cluster", seq_len(n_clades), "_ids")

  # ---- Column annotation (aligned to columns of M) ----
  # Build annotation in REVERSE order since pheatmap displays first column at bottom
  # Start with variant annotation (will be on bottom - first column)
  ann_col <- data.frame(Variant = dat$Variant, row.names = sample_ids)

  # Add T-stage annotation if requested (will be in middle - second column)
  if (annotate_t_stage && has_t) {
    # Turn T_stage into character - no need to handle NA since dataset is complete
    t_raw <- as.character(dat$T_stage)

    # Handle both individual T stages (T1, T2, T3, T4) and binned stages (T1-T2, T3-T4)
    # If values already start with "T" (e.g., "T3" or "T1-T2"), keep them;
    # if they look numeric (e.g., "3"), prepend "T".
    needs_T <- grepl("^[0-9]+$", t_raw)
    t_clean <- ifelse(needs_T, paste0("T", t_raw), t_raw)

    # Factor with canonical order (no Unknown needed since dataset is complete)
    # First check what format we have in the data
    unique_vals <- unique(t_clean)
    if (any(grepl("-", unique_vals))) {
      # We have binned values like "T1-T2", "T3-T4"
      t_lvls <- c("T1-T2", "T3-T4")
    } else {
      # We have individual values like "T1", "T2", "T3", "T4"
      t_lvls <- c("T1", "T2", "T3", "T4")
    }
    
    ann_col$`T Stage` <- factor(t_clean, levels = t_lvls)  # Use "T Stage" as column name for display
  }
  
  # Add cluster annotation last (will be on top - last column)
  # Create a mapping from sample_ids to final cluster assignments
  cluster_mapping <- setNames(cluster_df$Cluster, cluster_df$Patient_ID)
  original_ids <- unname(id_map[sample_ids])  # Convert sample_ids back to original Patient_IDs
  cluster_labels <- paste0("Cluster ", cluster_mapping[original_ids])
  names(cluster_labels) <- sample_ids
  
  # Use standard factor levels - legend order controlled by color order
  ann_col$Cluster <- factor(cluster_labels[sample_ids], levels = c("Cluster 1", "Cluster 2"))
  ann_col$Variant <- dat$Variant

  # Reorder rows of ann_col to match M's columns
  ann_col <- ann_col[colnames(M), , drop = FALSE]

  # ---- Annotation color lists ----
  # Build colors in same order as ann_col data frame (pheatmap displays in reverse)
  ann_colors <- list()
  
  # Add variant colors first (matches first column - will display at bottom)
  ann_colors$Variant <- variant_colors

  if (annotate_t_stage && has_t) {
    # Seed T-stage colors if not provided
    if (is.null(T_stage_colors)) T_stage_colors <- c()
    
    
    # Get the levels that actually appear in the data
    needed <- levels(ann_col$`T Stage`)  # Use "T Stage" column name
    
    # Find any missing colors and assign defaults
    missing <- setdiff(needed, names(T_stage_colors))
    if (length(missing)) {
      # Check if we have binned values and provide appropriate defaults
      if (any(grepl("-", missing))) {
        # Default colors for binned values
        default_binned <- c("T1-T2" = "#4575b4", "T3-T4" = "#d73027")
        fill_cols <- default_binned[missing]
        # For any still missing, use grays
        still_missing <- missing[is.na(fill_cols)]
        if (length(still_missing)) {
          gray_cols <- grDevices::gray.colors(length(still_missing), start = 0.3, end = 0.7)
          names(gray_cols) <- still_missing
          fill_cols[still_missing] <- gray_cols
        }
      } else {
        # Default for individual T stages
        fill_cols <- grDevices::gray.colors(length(missing), start = 0.3, end = 0.7)
        names(fill_cols) <- missing
      }
      T_stage_colors <- c(T_stage_colors, fill_cols)
    }
    
    ann_colors$`T Stage` <- T_stage_colors[needed]  # Use "T Stage" as key
  }
  
  # Add cluster colors last (matches last column - will display at top)
  # Order colors so Cluster 1 appears on top in legend
  if (!is.null(cluster_colors)) {
    # Keep the standard order of cluster colors
    ann_colors$Cluster <- cluster_colors[c("Cluster 1", "Cluster 2")]
  } else {
    # Default cluster colors if not provided
    default_cluster_colors <- c("Cluster 1" = "#94001E", "Cluster 2" = "#03507D")
    ann_colors$Cluster <- default_cluster_colors
  }
  
  # Add variant colors last (bottom annotation)
  ann_colors$Variant <- variant_colors

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
    group = group,
    feature_selector = feature_selector,
    top_features = top_features,
    ann_col = ann_col,
    ann_colors = ann_colors,
    cluster_df = cluster_df,  # Updated from clade_df
    clusters = clades,        # Updated from clades
    cluster_lists = cluster_lists,  # Updated from clade_lists
    heatmap_plot = heatmap_plot  # Plot object for patchwork
  )
}

# ===== File: Analysis/make_PCA.R =====
#' Create PCA plot with ellipses
#'
#' @param data Data frame with Patient_ID, Variant, and feature columns
#' @param plot_title Optional title for the plot (default: "")
#' @param ellipse_colors Named vector of colors for each variant/group
#' @param point_size Size of the points (default: 3 for standalone, 0.5 for multi-panel)
#' @param comp_x Which principal component to plot on x-axis (default: 1)
#' @param comp_y Which principal component to plot on y-axis (default: 2)
#' @return List containing the plot, PCA object, scores, scores_df, and explained variance
#' @export
make_PCA <- function(data, plot_title = "", 
                        ellipse_colors = c("PTC" = "#DF8D0A", "FV-PTC" = "#23744E", "FTC" = "#194992"),
                        point_size = 3, comp_x = 1, comp_y = 2) {
      #_Data preparation
      df <- as.data.frame(data)
      cls_col <- if ("Variant" %in% names(df)) "Variant" else names(df)[2]
      X <- df[, -c(1, 2), drop = FALSE]
      #_Coerce to numeric safely
      X[] <- lapply(X, function(v) suppressWarnings(as.numeric(v)))
      #_Handle NAs with median imputation
      if (anyNA(X)) {
        X[] <- lapply(X, function(v) {
          v[is.na(v)] <- stats::median(v, na.rm = TRUE)
          v
        })
      }
      Y <- factor(df[[cls_col]])
      
      #_Perform PCA
      pca <- stats::prcomp(X, center = TRUE, scale. = TRUE)
      
      #_Validate component indices
      max_comp <- min(ncol(X), nrow(X) - 1)
      if (comp_x > max_comp || comp_y > max_comp) {
        stop(paste("Requested components exceed available components. Max components:", max_comp))
      }
      
      scores <- pca$x[, c(comp_x, comp_y), drop = FALSE]
      explained <- round((pca$sdev^2 / sum(pca$sdev^2))[c(comp_x, comp_y)] * 100)
      
      #_Prepare plot data
      scores_df <- data.frame(
        Comp1 = scores[, 1],
        Comp2 = scores[, 2],
        Class = Y
      )
      
      # Identify NA values and create separate datasets
      na_mask <- is.na(scores_df$Class)
      scores_df_complete <- scores_df[!na_mask, , drop = FALSE]
      scores_df_na <- scores_df[na_mask, , drop = FALSE]
      
      #_Create PCA plot
      pca_plot <- ggplot2::ggplot() +
        # Plot complete cases with colors and ellipses
        {if(nrow(scores_df_complete) > 0) {
          list(
            ggplot2::geom_point(
              data = scores_df_complete,
              ggplot2::aes(x = Comp1, y = Comp2, color = Class, fill = Class),
              size = point_size, shape = 16
            ),
            ggplot2::stat_ellipse(
              data = scores_df_complete,
              ggplot2::aes(x = Comp1, y = Comp2, fill = Class),
              geom = "polygon", alpha = 0.3, color = NA
            )
          )
        }} +
        # Plot NA values as open circles without color
        {if(nrow(scores_df_na) > 0) {
          ggplot2::geom_point(
            data = scores_df_na,
            ggplot2::aes(x = Comp1, y = Comp2),
            size = point_size, shape = 1, color = "black", fill = NA
          )
        }} +
        ggplot2::scale_color_manual(values = ellipse_colors, drop = TRUE, na.translate = FALSE) +
        ggplot2::scale_fill_manual(values = ellipse_colors, drop = TRUE, na.translate = FALSE) +
        ggplot2::theme_minimal(base_family = "Arial") +
        ggplot2::labs(
          x = paste0("PC", comp_x, " (", explained[1], "%)"),
          y = paste0("PC", comp_y, " (", explained[2], "%)")
        ) +
        ggplot2::theme(
          axis.title = ggplot2::element_text(size = 25, face = "bold"),
          axis.text = ggplot2::element_text(size = 22, face = "bold", color = "black"),
          legend.position = "none",
          panel.grid.major = ggplot2::element_line(color = "gray80", linewidth = 0.8, linetype = "solid"),
          panel.grid.minor = ggplot2::element_blank(),
          panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 3.2),
          panel.background = ggplot2::element_blank()
        )
      
      #_Return useful objects for further analysis
      return(list(
        plot = pca_plot,
        pca = pca,
        scores = scores,
        scores_df = scores_df,
        explained = explained
      ))
}

# ===== File: Analysis/mummichog_ttests.R =====
#' Perform t-tests between groups and format results for Mummichog analysis
#'
#' @param data Data frame with Patient_ID, Variant, and feature columns
#' @param group_assignments Data frame with Patient_ID and group assignments (e.g., Clade)
#' @param group_column Name of the group column in group_assignments (default: "Clade")
#' @param output_filename Filename for the exported CSV results
#' @param group1_value Value representing the first group (default: 1)
#' @param group2_value Value representing the second group (default: 2)
#' @return List containing the results tibble and summary statistics
#' @export
mummichog_ttests <- function(data,
                            group_assignments,
                            group_column = "Clade",
                            output_filename,
                            group1_value = 1,
                            group2_value = 2) {
  
  #_Prepare data for t-tests with group assignments
  ttest_data <- data %>%
    dplyr::left_join(
      group_assignments %>% dplyr::mutate(Group_Test = !!rlang::sym(group_column)),
      by = "Patient_ID"
    ) %>%
    dplyr::select(-Patient_ID) %>%
    dplyr::filter(!is.na(Group_Test))
  
  #_Remove the original grouping column if it exists (and is different from Group_Test)
  if (group_column %in% names(ttest_data) && group_column != "Group_Test") {
    ttest_data <- ttest_data %>% dplyr::select(-!!rlang::sym(group_column))
  }
  
  #_Get feature names (exclude Group_Test and any remaining non-numeric columns)
  feature_names <- names(ttest_data)[names(ttest_data) != "Group_Test"]
  
  #_Filter to only include columns that look like metabolite features (HILIC or C18 prefix)
  metabolite_feature_names <- feature_names[stringr::str_starts(feature_names, "HILIC|C18")]
  
  #_Remove any remaining categorical columns that shouldn't be tested
  #_Only keep columns that can be converted to numeric for t-tests
  numeric_feature_names <- c()
  for (col in metabolite_feature_names) {
    test_values <- ttest_data[[col]]
    if (is.numeric(test_values) || (!is.factor(test_values) && !all(is.na(suppressWarnings(as.numeric(as.character(test_values))))))) {
      numeric_feature_names <- c(numeric_feature_names, col)
    }
  }
  feature_names <- numeric_feature_names
  
  cat("Total columns in data:", ncol(ttest_data), "\n")
  cat("Metabolite feature columns (HILIC/C18):", length(metabolite_feature_names), "\n")
  cat("Numeric feature columns for t-tests:", length(feature_names), "\n")
  
  #_Initialize results list
  ttest_results <- list()
  
  #_Perform t-test for each feature between groups
  cat("Performing t-tests for", length(feature_names), "features...\n")
  
  for (i in seq_along(feature_names)) {
    feature <- feature_names[i]
    
    #_Progress indicator
    if (i %% 5000 == 0) {
      cat("Processed", i, "of", length(feature_names), "features...\n")
    }
    
    group1_values <- ttest_data[ttest_data$Group_Test == group1_value, feature]
    group2_values <- ttest_data[ttest_data$Group_Test == group2_value, feature]
    
    #_Convert to numeric if factor and remove NA values
    if (is.factor(group1_values)) group1_values <- as.numeric(as.character(group1_values))
    if (is.factor(group2_values)) group2_values <- as.numeric(as.character(group2_values))
    group1_values <- group1_values[!is.na(group1_values)]
    group2_values <- group2_values[!is.na(group2_values)]
    
    #_Check if there's sufficient data for t-test
    if (length(group1_values) < 2 || length(group2_values) < 2) {
      ttest_results[[feature]] <- 1.0  # Assign p-value of 1 for insufficient data
      next
    }
    
    #_Check for constant data (no variance) - assign p-value of 1
    #_Use a safer method to check for constant data
    group1_constant <- length(unique(group1_values)) == 1
    group2_constant <- length(unique(group2_values)) == 1
    if (group1_constant && group2_constant) {
      ttest_results[[feature]] <- 1.0  # Assign p-value of 1 for constant data
      next
    }
    
    #_Perform t-test with error handling
    tryCatch({
      test_result <- t.test(group1_values, group2_values)
      ttest_results[[feature]] <- test_result$p.value
    }, error = function(e) {
      cat("Warning: Assigning p-value of 1 to feature", feature, "due to error:", e$message, "\n")
      ttest_results[[feature]] <- 1.0  # Assign p-value of 1 for errors
    })
  }
  
  #_Parse feature names and create results tibble
  results_tibble <- tibble::tibble(
    Feature = names(ttest_results),
    p.value = unlist(ttest_results)
  ) %>%
    dplyr::mutate(
      #_Extract mode and convert: HILIC -> pos, C18 -> neg
      mode = dplyr::case_when(
        stringr::str_starts(Feature, "HILIC") ~ "positive",
        stringr::str_starts(Feature, "C18") ~ "negative",
        TRUE ~ NA_character_
      ),
      #_Extract m.z (first number after underscore)
      m.z = stringr::str_extract(Feature, "(?<=_)[0-9.]+"),
      #_Extract rt (second number - after second underscore)
      r.t = stringr::str_extract(Feature, "_[0-9.]+_([0-9.]+)") %>%
        stringr::str_extract("[0-9.]+$")
    ) %>%
    dplyr::select(m.z, p.value, mode, r.t) %>%
    dplyr::mutate(
      m.z = as.numeric(m.z),
      r.t = as.numeric(r.t)
    ) %>%
    #_Remove rows where feature parsing failed (invalid feature names)
    dplyr::filter(!is.na(m.z) & !is.na(mode))
  
  #_Export results
  readr::write_csv(results_tibble, paste0("Outputs/Mummichog Inputs/", output_filename))
  
  #_Display summary
  cat("T-test results exported to:", paste0("Outputs/", output_filename), "\n")
  cat("Total features processed:", length(feature_names), "\n")
  cat("Features with actual t-test p-values:", sum(unlist(ttest_results) < 1.0), "\n")
  cat("Features with assigned p-value of 1 (constant/insufficient data):", sum(unlist(ttest_results) == 1.0), "\n")
  cat("Total results:", nrow(results_tibble), "\n")
  cat("Group", group1_value, "vs Group", group2_value, "comparison\n")
  cat("First 10 results:\n")
  print(head(results_tibble, 10))
  
  #_Return results
  return(list(
    results = results_tibble,
    n_features = nrow(results_tibble),
    output_file = paste0("Outputs/", output_filename),
    group1_value = group1_value,
    group2_value = group2_value
  ))
}

# ===== File: Analysis/targeted_comp.R =====
#' Perform statistical comparisons on targeted metabolomics data
#'
#' @param data Data frame with metabolite data, must include Patient_ID and grouping variable
#' @param grouping_var Character string specifying the column name for grouping (e.g., "Clade", "Variant", "T")
#' @param test_type Character string specifying test type: "t_test" or "anova"
#' @param exclude_cols Character vector of column names to exclude from analysis (default: c("Patient_ID"))
#' @param exclude_isotopes Logical, whether to exclude isotope metabolites (15N, 13C) from analysis (default: TRUE)
#' @return Tibble with statistical test results including p-values, FDR correction, and group means
#' @export
targeted_metabolite_comparison <- function(data, grouping_var, test_type = "t_test", exclude_cols = c("Patient_ID"), exclude_isotopes = TRUE) {
  
  # Load required libraries
  library(dplyr)
  library(purrr)
  library(broom)
  
  # Identify metabolite columns
  metabolite_cols <- names(data)[!names(data) %in% c(exclude_cols, grouping_var)]
  
  # Filter out isotopes if requested
  if (exclude_isotopes) {
    metabolite_cols <- metabolite_cols[!str_detect(metabolite_cols, "15N|13C")]
  }
  
  # Check if grouping variable exists
  if (!grouping_var %in% names(data)) {
    stop(paste("Grouping variable", grouping_var, "not found in data"))
  }
  
  # Perform statistical tests for each metabolite
  test_results <- map_dfr(metabolite_cols, ~{
    # Extract data for current metabolite
    metabolite_data <- data %>%
      filter(!is.na(.data[[grouping_var]]), !is.na(.data[[.x]])) %>%
      select(group = all_of(grouping_var), metabolite = all_of(.x))
    
    # Convert group to factor
    metabolite_data$group <- as.factor(metabolite_data$group)
    
    # Check if we have sufficient data
    if(nrow(metabolite_data) > 0 && 
       length(unique(metabolite_data$group)) >= 2 &&
       all(table(metabolite_data$group) >= 2)) {
      
      if(test_type == "t_test") {
        # Perform t-test (only works for 2 groups)
        if(length(unique(metabolite_data$group)) == 2) {
          test_result <- t.test(metabolite ~ group, data = metabolite_data)
          
          # Calculate group means
          group_means <- metabolite_data %>%
            group_by(group) %>%
            summarise(mean_val = mean(metabolite, na.rm = TRUE), .groups = "drop")
          
          # Create named list of means
          means_list <- setNames(group_means$mean_val, paste0("mean_", group_means$group))
          
          # Calculate fold change (difference for log2 data)
          if(nrow(group_means) == 2) {
            fold_change <- group_means$mean_val[2] - group_means$mean_val[1]
          } else {
            fold_change <- NA
          }
          
          result <- tibble(
            Metabolite = .x,
            test_statistic = test_result$statistic,
            p_value = test_result$p.value,
            fold_change = fold_change
          )
          
          # Add group means as separate columns
          for(i in seq_along(means_list)) {
            result[[names(means_list)[i]]] <- means_list[[i]]
          }
          
        } else {
          # More than 2 groups - cannot perform t-test
          result <- tibble(
            Metabolite = .x,
            test_statistic = NA,
            p_value = NA,
            fold_change = NA,
            error = "t-test requires exactly 2 groups"
          )
        }
        
      } else if(test_type == "anova") {
        # Perform ANOVA
        test_result <- aov(metabolite ~ group, data = metabolite_data)
        anova_summary <- summary(test_result)
        
        # Calculate group means
        group_means <- metabolite_data %>%
          group_by(group) %>%
          summarise(mean_val = mean(metabolite, na.rm = TRUE), .groups = "drop")
        
        # Create named list of means
        means_list <- setNames(group_means$mean_val, paste0("mean_", group_means$group))
        
        # Calculate overall range (max - min means)
        fold_change <- max(group_means$mean_val) - min(group_means$mean_val)
        
        result <- tibble(
          Metabolite = .x,
          test_statistic = anova_summary[[1]]$`F value`[1],
          p_value = anova_summary[[1]]$`Pr(>F)`[1],
          fold_change = fold_change
        )
        
        # Add group means as separate columns
        for(i in seq_along(means_list)) {
          result[[names(means_list)[i]]] <- means_list[[i]]
        }
        
      } else {
        stop("test_type must be either 't_test' or 'anova'")
      }
      
    } else {
      # Insufficient data
      result <- tibble(
        Metabolite = .x,
        test_statistic = NA,
        p_value = NA,
        fold_change = NA,
        error = "Insufficient data"
      )
    }
    
    return(result)
  })
  
  # Add FDR correction and significance flags
  final_results <- test_results %>%
    mutate(
      p_value_fdr = p.adjust(p_value, method = "fdr"),
      significant = p_value < 0.05 & !is.na(p_value),
      significant_fdr = p_value_fdr < 0.05 & !is.na(p_value_fdr)
    ) %>%
    arrange(p_value_fdr)
  
  # Report results
  n_significant <- sum(final_results$significant, na.rm = TRUE)
  n_significant_fdr <- sum(final_results$significant_fdr, na.rm = TRUE)
  
  cat("Statistical comparison results:\n")
  cat("Test type:", test_type, "\n")
  cat("Grouping variable:", grouping_var, "\n")
  cat("Isotopes excluded:", exclude_isotopes, "\n")
  cat("Number of metabolites tested:", nrow(final_results), "\n")
  cat("Significantly different (uncorrected p < 0.05):", n_significant, "\n")
  cat("Significantly different (FDR corrected p < 0.05):", n_significant_fdr, "\n\n")
  
  return(final_results)
}

# ===== File: Other/combine_R_files.r =====
combine_R_files <- function(input_dir, output_file) {
  input_dir <- normalizePath(input_dir, mustWork = TRUE)

  # Find all .R or .r files recursively (case-insensitive)
  r_files <- list.files(
    input_dir,
    pattern = "(?i)\\.r$", full.names = TRUE, recursive = TRUE
  )

  if (length(r_files) == 0) {
    stop("❌ No .R or .r files found in ", input_dir)
  }

  # Natural sort (00, 01, 02 … instead of 1, 10, 2)
  if (!requireNamespace("gtools", quietly = TRUE)) {
    install.packages("gtools")
  }
  r_files <- gtools::mixedsort(r_files)

  # Get relative paths for cleaner headers
  rel_paths <- sub(paste0("^", input_dir, "/?"), "", r_files)

  # Read and combine with headers
  all_contents <- unlist(Map(function(f, rel) {
    c(
      paste0("# ===== File: ", rel, " ====="),
      readLines(f, warn = FALSE),
      ""
    )
  }, r_files, rel_paths))

  # Count total number of lines
  total_lines <- length(all_contents)

  # Ensure output directory exists
  out_dir <- dirname(output_file)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  # Write out combined file
  writeLines(all_contents, output_file)
  message(
    "✅ Combined ", length(r_files), " R files into: ", normalizePath(output_file),
    "\n📏 Total lines in combined file: ", total_lines
  )
}

# ===== File: Other/list_tree.r =====
list_tree <- function(path = ".", prefix = "") {
  items <- list.files(path, full.names = TRUE)
  for (i in seq_along(items)) {
    item <- items[i]
    name <- basename(item)
    cat(prefix, if (i == length(items)) "└── " else "├── ", name, "\n", sep = "")
    if (dir.exists(item)) {
      list_tree(item, paste0(prefix, if (i == length(items)) "    " else "│   "))
    }
  }
}

# ===== File: Other/load_checkpoints.r =====
#' Load all checkpoint .rds files under raw_path/Checkpoints
#' @param raw_path Path containing the "Checkpoints" directory
#' @param envir Environment to load into (default = .GlobalEnv)
load_checkpoints <- function(raw_path, envir = .GlobalEnv) {
  checkpoint_dir <- file.path(raw_path, "Checkpoints")

  if (!dir.exists(checkpoint_dir)) {
    stop("❌ Checkpoint directory not found: ", checkpoint_dir)
  }

  # find all .rds files recursively
  rds_files <- list.files(
    checkpoint_dir,
    pattern = "\\.rds$",
    full.names = TRUE,
    recursive = TRUE
  )

  if (length(rds_files) == 0) {
    stop("❌ No .rds files found in ", checkpoint_dir)
  }

  loaded_names <- c()

  for (f in rds_files) {
    message(">>> Loading checkpoint: ", f)
    objs <- readRDS(f)
    list2env(objs, envir = envir)
    loaded_names <- c(loaded_names, names(objs))
  }

  loaded_names <- unique(loaded_names)
  message(
    "✅ Loaded ", length(loaded_names), " unique objects from ",
    length(rds_files), " checkpoint file(s)"
  )

  invisible(loaded_names)
}

# ===== File: Other/load_with_checkpoints.r =====
load_with_checkpoints <- function(raw_path, files_csv, files_xlsx, checkpoint_dir = here::here("R", "Checkpoints")) {
  # Ensure checkpoint directory exists
  if (!dir.exists(checkpoint_dir)) dir.create(checkpoint_dir, recursive = TRUE)
  rds_file <- file.path(checkpoint_dir, "all_objects.rds")

  if (file.exists(rds_file)) {
    message(">>> Loading from checkpoint: ", rds_file)
    objs <- readRDS(rds_file)
    list2env(objs, envir = .GlobalEnv)
    rm(objs)
  } else {
    message(">>> No checkpoint found, reading raw data...")

    # --- Load CSVs ---
    csv_objs <- purrr::map2(
      files_csv$name, files_csv$file,
      ~ readr::read_csv(file.path(raw_path, paste0(.y, ".csv")))
    )
    names(csv_objs) <- files_csv$name

    # --- Load XLSX ---
    xlsx_objs <- purrr::pmap(
      list(files_xlsx$name, files_xlsx$file, files_xlsx$sheet),
      ~ readxl::read_excel(file.path(raw_path, paste0(..2, ".xlsx")), sheet = ..3)
    )
    names(xlsx_objs) <- files_xlsx$name

    # Combine all objects
    all_objs <- c(csv_objs, xlsx_objs)

    # Save checkpoint
    saveRDS(all_objs, rds_file)
    message(">>> Saved temporary checkpoint: ", rds_file)

    # Expose to environment
    list2env(all_objs, envir = .GlobalEnv)
  }
}

# ===== File: Other/promote_checkpoint.r =====
promote_checkpoint <- function(raw_path) {
  checkpoint_dir <- here::here("R", "Checkpoints")
  temp_rds <- file.path(checkpoint_dir, "all_objects.rds")
  final_rds <- file.path(raw_path, "all_objects.rds")

  if (file.exists(temp_rds)) {
    file.copy(temp_rds, final_rds, overwrite = TRUE)
    message(">>> Final checkpoint saved to: ", final_rds)
    unlink(checkpoint_dir, recursive = TRUE, force = TRUE)
    message(">>> Temporary checkpoint directory deleted: ", checkpoint_dir)
  } else {
    message("⚠️ No temporary checkpoint found at end of pipeline")
  }
}

# ===== File: Other/run_scripts.r =====
run_scripts <- function(spec = "06", raw_path) {
  purrr::walk(
    list.files(
      here::here("R", "Utilities"),
      pattern = "\\.[rR]$",
      full.names = TRUE,
      recursive = TRUE
    ),
    source
  )

  script_dir <- here::here("R", "Scripts")
  scripts <- gtools::mixedsort(list.files(script_dir, pattern = "\\.[rR]$", full.names = TRUE))

  # Put checkpoints inside raw_path
  chk_dir <- file.path(raw_path, "Checkpoints")
  if (!dir.exists(chk_dir)) dir.create(chk_dir, recursive = TRUE)

  # Helper: extract prefix (e.g. "01", "02")
  script_ids <- substr(basename(scripts), 1, 2)

  # Parse spec argument
  to_run <- c()
  parts <- unlist(strsplit(spec, ","))
  for (p in parts) {
    if (grepl("-", p)) {
      range <- unlist(strsplit(p, "-"))
      start <- which(script_ids == range[1])
      end <- which(script_ids == range[2])
      if (length(start) == 0 | length(end) == 0) stop("Invalid range: ", p)
      to_run <- c(to_run, seq(start, end))
    } else {
      idx <- which(script_ids == p)
      if (length(idx) == 0) stop("Script not found for ID: ", p)
      to_run <- c(to_run, idx)
    }
  }
  to_run <- unique(sort(to_run))

  # Run scripts
  for (i in to_run) {
    script <- scripts[i]
    script_name <- tools::file_path_sans_ext(basename(script))
    message(">>> Running: ", basename(script))

    tryCatch(
      {
        source(script, local = .GlobalEnv, echo = TRUE, keep.source = TRUE)

        # Save checkpoint after each script
        chk_file <- file.path(chk_dir, paste0(script_name, ".rds"))
        saveRDS(as.list(.GlobalEnv), chk_file)
        message(">>> Checkpoint saved: ", chk_file)
      },
      error = function(e) {
        message("❌ Error in ", basename(script), ": ", e$message)
        stop(e)
      }
    )
  }

  message("✅ Finished running scripts: ", spec)
}

# ===== File: Preprocessing/assign_T_stage.R =====
#' Assign AJCC 8th-edition T stage (T1–T4) for differentiated thyroid cancer
#'
#' Determines the primary tumor T category using largest tumor dimension (LD)
#' and extrathyroidal extension (ETE) per AJCC/UICC TNM 8th edition.
#'
#' @section Assumption:
#' ETE is coded as 0 = none, 1 = microscopic/minimal, 2 = gross/extensive.
#' This function assumes ETE = 2 always represents gross ETE beyond strap muscles
#' and therefore upgrades to at least T4 if size-based T-stage is lower.
#'
#' @param df Data frame containing tumor size and ETE columns.
#' @param ld_col Character. Column name with largest tumor dimension in cm (default "LD").
#' @param ete_col Character. Column name with ETE code 0/1/2 (default "ETE").
#' @param units Either "cm" (default) or "mm" for the size column.
#' @param out_col Character. Name of output column to write (default "T_stage").
#'
#' @details
#' Size-based rules (microscopic ETE ignored in AJCC8):
#' \itemize{
#'   \item T1: \eqn{\le} 2 cm, confined to thyroid
#'   \item T2: > 2 cm and \eqn{\le} 4 cm, confined to thyroid
#'   \item T3: > 4 cm, confined to thyroid
#' }
#' ETE handling:
#' \itemize{
#'   \item ETE == 0 or 1: keep size-based T (microscopic ETE does not upstage)
#'   \item ETE == 2: upgrade to T4 if size-based T-stage is lower
#' }
#' Lymphovascular invasion (LVI) and multifocality are not used for T1–T4 assignment.
#'
#' @return A copy of \code{df} with an ordered factor column \code{out_col} in \{T1,T2,T3,T4\}.
#' @export
assign_T_stage <- function(
    df,
    ld_col = "LD",
    ete_col = "ETE",
    units = c("cm", "mm"),
    out_col = "T_stage") {
  units <- match.arg(units)

  LD <- df[[ld_col]]
  ETE <- df[[ete_col]]

  # Convert to cm if needed
  LD_cm <- if (units == "mm") LD / 10 else LD

  # Size-only T (AJCC8; microscopic ETE ignored)
  T_size <- dplyr::case_when(
    is.na(LD_cm) ~ NA_character_,
    LD_cm <= 2 ~ "T1",
    LD_cm > 2 & LD_cm <= 4 ~ "T2",
    LD_cm > 4 ~ "T3",
    TRUE ~ NA_character_
  )

  # Upgrade logic for gross ETE (ETE==2 → always T4 if lower)
  ord <- c(T1 = 1L, T2 = 2L, T3 = 3L, T4 = 4L)
  T_final <- dplyr::case_when(
    is.na(T_size) ~ NA_character_,
    is.na(ETE) ~ T_size,
    ETE %in% c(0, 1) ~ T_size,
    ETE == 2 ~ {
      cur_ord <- ord[T_size]
      ifelse(cur_ord < ord["T4"], "T4", T_size)
    },
    TRUE ~ T_size
  )

  df[[out_col]] <- factor(T_final, levels = c("T1", "T2", "T3", "T4"), ordered = TRUE)
  df
}

# ===== File: Preprocessing/clean_pathway_names.R =====
clean_pathway_names <- function(pathway_names) {
  pathway_names %>%
    # Handle Greek letters and special cases first
    stringr::str_replace_all(stringr::regex("\\bBeta[- ]Alanine\\b", ignore_case = TRUE), "β-Alanine") %>%
    stringr::str_replace_all(stringr::regex("\\balpha[- ]?linolenic\\b", ignore_case = TRUE), "α-Linolenic") %>%
    stringr::str_replace_all(stringr::regex("\\bgamma[- ]?linol", ignore_case = TRUE), "γ-linol") %>%
    stringr::str_replace_all(stringr::regex("\\bgama[- ]?linoleic\\b", ignore_case = TRUE), "γ-Linoleic") %>%
    # Vitamin formatting - keep ONLY parenthetical content
    stringr::str_replace_all(stringr::regex("\\bVitamin A \\(Retinol\\)\\b", ignore_case = TRUE), "Retinol") %>%
    stringr::str_replace_all(stringr::regex("\\bVitamin B1 \\(thiamin\\)\\b", ignore_case = TRUE), "Thiamin") %>%
    stringr::str_replace_all(stringr::regex("\\bVitamin D3 \\(cholecalciferol\\)\\b", ignore_case = TRUE), "Cholecalciferol") %>%
    # Amino acid abbreviations (excluding when part of beta-alanine or other compounds)
    stringr::str_replace_all(stringr::regex("\\bAlanine(?!\\s)", ignore_case = TRUE), "Ala") %>%
    stringr::str_replace_all(stringr::regex("(?<!β-)\\bAlanine\\b", ignore_case = TRUE), "Ala") %>%
    stringr::str_replace_all(stringr::regex("\\bArginine\\b", ignore_case = TRUE), "Arg") %>%
    stringr::str_replace_all(stringr::regex("\\bAsparagine\\b", ignore_case = TRUE), "Asn") %>%
    stringr::str_replace_all(stringr::regex("\\bAspartate\\b", ignore_case = TRUE), "Asp") %>%
    stringr::str_replace_all(stringr::regex("\\bCysteine\\b", ignore_case = TRUE), "Cys") %>%
    stringr::str_replace_all(stringr::regex("\\bGlutamate\\b", ignore_case = TRUE), "Glu") %>%
    stringr::str_replace_all(stringr::regex("\\bGlutamine\\b", ignore_case = TRUE), "Gln") %>%
    stringr::str_replace_all(stringr::regex("\\bGlycine\\b", ignore_case = TRUE), "Gly") %>%
    stringr::str_replace_all(stringr::regex("\\bHistidine\\b", ignore_case = TRUE), "His") %>%
    stringr::str_replace_all(stringr::regex("\\bIsoleucine\\b", ignore_case = TRUE), "Ile") %>%
    stringr::str_replace_all(stringr::regex("\\bLeucine\\b", ignore_case = TRUE), "Leu") %>%
    stringr::str_replace_all(stringr::regex("\\bLysine\\b", ignore_case = TRUE), "Lys") %>%
    stringr::str_replace_all(stringr::regex("\\bMethionine\\b", ignore_case = TRUE), "Met") %>%
    stringr::str_replace_all(stringr::regex("\\bPhenylalanine\\b", ignore_case = TRUE), "Phe") %>%
    stringr::str_replace_all(stringr::regex("\\bProline\\b", ignore_case = TRUE), "Pro") %>%
    stringr::str_replace_all(stringr::regex("\\bSerine\\b", ignore_case = TRUE), "Ser") %>%
    stringr::str_replace_all(stringr::regex("\\bThreonine\\b", ignore_case = TRUE), "Thr") %>%
    stringr::str_replace_all(stringr::regex("\\bTryptophan\\b", ignore_case = TRUE), "Trp") %>%
    stringr::str_replace_all(stringr::regex("\\bTyrosine\\b", ignore_case = TRUE), "Tyr") %>%
    stringr::str_replace_all(stringr::regex("\\bValine\\b", ignore_case = TRUE), "Val") %>%
    stringr::str_replace_all(stringr::regex("\\bTaurine\\b", ignore_case = TRUE), "Tau") %>%
    stringr::str_replace_all(stringr::regex("\\bHypotaurine\\b", ignore_case = TRUE), "Hypotau") %>%
    # Fatty acid abbreviations
    stringr::str_replace_all(stringr::regex("\\bunsaturated fatty acids?\\b", ignore_case = TRUE), "UFAs") %>%
    # Prostaglandin abbreviation
    stringr::str_replace_all(stringr::regex("\\bProstaglandin\\b", ignore_case = TRUE), "PG") %>%
    # Special abbreviations and replacements
    stringr::str_replace_all(stringr::regex("\\bcytochrome P450\\b", ignore_case = TRUE), "Cyp450") %>%
    stringr::str_replace_all(stringr::regex("\\bGlycosylphosphatidylinositol \\(GPI\\)-anchor\\b", ignore_case = TRUE), "GPI-anchor") %>%
    stringr::str_replace_all(stringr::regex("\\bThiamine\\b", ignore_case = TRUE), "Vitamin B1") %>%
    # Replace 'and' with '&'
    stringr::str_replace_all(stringr::regex("\\band\\b", ignore_case = FALSE), " & ") %>%
    # Universal capitalizations - specific patterns first
    stringr::str_replace_all(stringr::regex("\\bbiosynthesis\\s*&\\s*metabolism\\b", ignore_case = TRUE), "Biosynthesis/Metabolism") %>%
    stringr::str_replace_all(stringr::regex("\\bde\\s+novo\\s+fatty\\b", ignore_case = TRUE), "De Novo Fatty") %>%
    stringr::str_replace_all(stringr::regex("\\bperoxisomal\\s+oxidation\\b", ignore_case = TRUE), "Peroxisomal Oxidation") %>%
    stringr::str_replace_all(stringr::regex("\\bcarnitine\\s+shuttle\\b", ignore_case = TRUE), "Carnitine Shuttle") %>%
    stringr::str_replace_all(stringr::regex("\\bmetabolism\\b", ignore_case = TRUE), "Metabolism") %>%
    stringr::str_replace_all(stringr::regex("\\bdegradation\\b", ignore_case = TRUE), "Degradation") %>%
    stringr::str_replace_all(stringr::regex("\\bbiosynthesis\\b", ignore_case = TRUE), "Biosynthesis") %>%
    stringr::str_replace_all(stringr::regex("\\bshuttle\\b", ignore_case = TRUE), "Shuttle") %>%
    stringr::str_replace_all(stringr::regex("\\bnovo\\b", ignore_case = TRUE), "Novo") %>%
    stringr::str_replace_all(stringr::regex("\\bfatty\\b", ignore_case = TRUE), "Fatty") %>%
    stringr::str_replace_all(stringr::regex("\\bperoxisomal\\b", ignore_case = TRUE), "Peroxisomal") %>%
    stringr::str_replace_all(stringr::regex("\\bformation\\b", ignore_case = TRUE), "Formation") %>%
    stringr::str_replace_all(stringr::regex("\\barachidonate\\b", ignore_case = TRUE), "Arachidonate") %>%
    stringr::str_replace_all(stringr::regex("\\bdihomo\\b", ignore_case = TRUE), "Dihomo") %>%
    stringr::str_replace_all(stringr::regex("\\bactivation\\b", ignore_case = TRUE), "Activation") %>%
    stringr::str_replace_all(stringr::regex("\\bmannose\\b", ignore_case = TRUE), "Mannose") %>%
    stringr::str_replace_all(stringr::regex("\\bacid\\b", ignore_case = TRUE), "Acid") %>%
    stringr::str_replace_all(stringr::regex("\\boxidation\\b", ignore_case = TRUE), "Oxidation") %>%
    stringr::str_replace_all(stringr::regex("\\bperoxisome\\b", ignore_case = TRUE), "Peroxisome") %>%
    stringr::str_replace_all(stringr::regex("\\bretinol\\b", ignore_case = TRUE), "Retinol") %>%
    stringr::str_replace_all(stringr::regex("\\bbile\\b", ignore_case = TRUE), "Bile") %>%
    stringr::str_replace_all(stringr::regex("\\bhormone\\b", ignore_case = TRUE), "Hormone") %>%
    stringr::str_replace_all(stringr::regex("\\bother enzymes\\b", ignore_case = TRUE), "Other Enzymes") %>%
    stringr::str_replace_all(stringr::regex("\\bnicotinamide\\b", ignore_case = TRUE), "Nicotinamide") %>%
    stringr::str_replace_all(stringr::regex("\\bnicotinate\\b", ignore_case = TRUE), "Nicotinate") %>%
    stringr::str_replace_all(stringr::regex("\\bcycle\\b", ignore_case = TRUE), "Cycle") %>%
    stringr::str_replace_all(stringr::regex("\\bfolate\\b", ignore_case = TRUE), "Folate") %>%
    stringr::str_replace_all(stringr::regex("\\bsteroid\\b", ignore_case = TRUE), "Steroid") %>%
    stringr::str_replace_all(stringr::regex("\\bprimary\\b", ignore_case = TRUE), "Primary") %>%
    stringr::str_replace_all(stringr::regex("\\bsphingolipid\\b", ignore_case = TRUE), "Sphingolipid") %>%
    stringr::str_replace_all(stringr::regex("\\bporphyrin\\b", ignore_case = TRUE), "Porphyrin") %>%
    stringr::str_replace_all(stringr::regex("\\bglucuronate\\b", ignore_case = TRUE), "Glucuronate") %>%
    stringr::str_replace_all(stringr::regex("\\binterconversions\\b", ignore_case = TRUE), "Interconversions") %>%
    stringr::str_replace_all(stringr::regex("\\bascorbate\\b", ignore_case = TRUE), "Ascorbate") %>%
    stringr::str_replace_all(stringr::regex("\\baldarate\\b", ignore_case = TRUE), "Aldarate") %>%
    stringr::str_replace_all(stringr::regex("\\bdicarboxylate\\b", ignore_case = TRUE), "Dicarboxylate")
}

# ===== File: Preprocessing/process_feature_table.R =====
process_feature_table <- function(data_tibble, sequence_tibble, data_type) {
  data_tibble %>%
    mutate(Feature = paste(data_type, mz, time, sep = "_")) %>%
    select(Feature, everything(), -c(mz, time)) %>%
    rename_with(~ gsub("\\.mzXML$", "", .)) %>%
    pivot_longer(cols = -Feature, names_to = "File_Name", values_to = "Value") %>%
    pivot_wider(names_from = Feature, values_from = Value) %>%
    left_join(sequence_tibble, by = "File_Name") %>%
    select(Sample_ID, everything(), -File_Name) %>%
    mutate(Sample_ID = sub("_.*", "", Sample_ID)) %>%
    arrange(Sample_ID) %>%
    select(-Batch)
}

# ===== File: Preprocessing/t_and_name_FT.R =====
t_and_name_FT <- function(data_tibble, sequence_tibble) {
  data_tibble %>%
    pivot_longer(cols = -Feature, names_to = "File_Name", values_to = "Value") %>%
    pivot_wider(names_from = Feature, values_from = Value) %>%
    left_join(sequence_tibble, by = "File_Name") %>% # Join based on File_Name
    select(Sample_ID, everything(), -File_Name) %>% # Make Sample_ID the first column and remove File_Name
    arrange(Sample_ID)
}

# ===== File: Visualization/print_to_png.R =====
#' Print plot to PNG with auto-refresh for macOS Preview
#'
#' @param plot The plot object to print
#' @param filename Name of the PNG file (with or without .png extension)
#' @param width Width in inches (default: 8.5)
#' @param height Height in inches (default: 11)
#' @param dpi Resolution in DPI (default: 300 for high quality)
#' @param output_dir Directory to save the PNG (default: "Outputs")
#' @param auto_open Whether to automatically open in Preview on first run (default: TRUE)
#' @return Invisible path to the created PNG file
#' @export
print_to_png <- function(plot, filename, width = 8.5, height = 11, dpi = 600,
                         output_dir = "Figures", auto_open = TRUE) {
  # Ensure filename has .png extension
  if (!grepl("\\.png$", filename, ignore.case = TRUE)) {
    filename <- paste0(filename, ".png")
  }

  # Create full path
  filepath <- file.path(output_dir, filename)

  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Check if file already exists (for auto-open logic)
  file_exists <- file.exists(filepath)

  # Save the plot as PNG
  ggplot2::ggsave(
    filename = filepath,
    plot = plot,
    width = width,
    height = height,
    dpi = dpi,
    units = "in",
    device = "png"
  )

  # Auto-open in Preview only on first run (or if specified)
  if (auto_open && !file_exists) {
    system(paste("open", shQuote(filepath)))
    cat("PNG saved and opened in Preview:", filepath, "\n")
    cat("Preview will auto-refresh when you re-run this function!\n")
  } else {
    cat("PNG updated:", filepath, "\n")
  }

  # Return path invisibly
  invisible(filepath)
}

# ===== File: Visualization/publication_theme.R =====
#* Base theme (never hides legend)
  theme_pub_base <- function(base_size = 12, base_family = "Arial",
                             border_linewidth = 0.8,
                             major_frac = 0.5, # major grid = 50% of border
                             minor_frac = 0.5, # minor grid = 50% of major (i.e., 25% of border)
                             show_major_x = TRUE, show_major_y = FALSE,
                             show_minor_x = FALSE, show_minor_y = FALSE) {
    major_lwd <- border_linewidth * major_frac
    minor_lwd <- major_lwd * minor_frac

    ggprism::theme_prism(base_size = base_size, base_family = base_family) +
      ggplot2::theme(
        plot.margin = grid::unit(c(8, 8, 6, 6), "pt"),
        axis.line = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_line(linewidth = 0.3),
        axis.ticks.length = grid::unit(0.2, "cm"),
        axis.title.x = ggplot2::element_text(size = 9, face = "bold"),
        axis.title.y = ggplot2::element_text(size = 9, face = "bold"),
        axis.text.x = ggplot2::element_text(size = 8, face = "bold"),
        axis.text.y = ggplot2::element_text(size = 8, face = "bold"),

        # Border
        panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = border_linewidth),
        panel.background = ggplot2::element_blank(),

        # Major grids
        panel.grid.major.x = if (isTRUE(show_major_x)) {
          ggplot2::element_line(color = "gray80", linewidth = major_lwd)
        } else {
          ggplot2::element_blank()
        },
        panel.grid.major.y = if (isTRUE(show_major_y)) {
          ggplot2::element_line(color = "gray80", linewidth = major_lwd)
        } else {
          ggplot2::element_blank()
        },

        # Minor grids
        panel.grid.minor.x = if (isTRUE(show_minor_x)) {
          ggplot2::element_line(color = "gray90", linewidth = minor_lwd)
        } else {
          ggplot2::element_blank()
        },
        panel.grid.minor.y = if (isTRUE(show_minor_y)) {
          ggplot2::element_line(color = "gray90", linewidth = minor_lwd)
        } else {
          ggplot2::element_blank()
        },
        legend.title = ggplot2::element_text(size = 9, face = "bold"),
        legend.text = ggplot2::element_text(size = 8, face = "bold")
      )
  }
#* Simple (no legend)
  theme_pub_simple <- function(base_size = 12, base_family = "Arial",
                               border_linewidth = 0.8,
                               major_frac = 0.5, minor_frac = 0.5) {
    theme_pub_base(base_size, base_family,
      border_linewidth = border_linewidth,
      major_frac = major_frac, minor_frac = minor_frac,
      show_major_x = TRUE, show_major_y = FALSE,
      show_minor_x = FALSE, show_minor_y = FALSE
    )
  }
#* Dot Bar
  theme_pub_dotbar <- function(base_size = 12, base_family = "Arial",
                              border_linewidth = 0.8,
                              bar_border_scale = 1.5) {
    theme_pub_base(
      base_size, base_family,
      border_linewidth = border_linewidth * bar_border_scale,
      major_frac = 0.25, minor_frac = 0.25,
      show_major_x = FALSE, show_major_y = FALSE,
      show_minor_x = FALSE, show_minor_y = FALSE # ⟵ only horizontal minor grid
    ) +
      ggplot2::theme(
        panel.ontop       = FALSE, # gridlines behind bars/dots
        axis.text.x       = ggplot2::element_text(size = base_size * 0.75, angle = 0, vjust = 1),
        legend.position   = "top",
        legend.direction  = "horizontal",
        legend.box        = "horizontal",
        legend.margin     = ggplot2::margin(t = 0, r = 0, b = -4, l = 0),
        legend.key.width  = grid::unit(0.25, "cm"),
        legend.key.height = grid::unit(0.25, "cm"),
        legend.key.size   = grid::unit(0.25, "cm"),
        legend.spacing.x  = grid::unit(0.05, "cm"),
        legend.title      = ggplot2::element_blank()
      )
  }
#* Grouped (keeps legend; choose position) ----
  theme_pub_grouped <- function(base_size = 12, base_family = "Arial",
                                legend_position = "right",
                                legend_direction = "vertical") {
    theme_pub_base(base_size, base_family) +
      ggplot2::theme(
        legend.position = legend_position,
        legend.direction = legend_direction,
        # match bar plot frame explicitly
        panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.5),
        panel.grid.major.x = ggplot2::element_line(color = "gray80", linewidth = 0.3),
        panel.grid.major.y = ggplot2::element_blank()
      )
  }
#* Bar plot
  theme_pub_barplot <- function() {
    theme_pub_dotbar() +
      theme(
        panel.ontop        = FALSE,
        # kill vertical gridlines:
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        # keep only horizontal major grid:
        panel.grid.major.y = element_line(color = "gray80", linewidth = 0.3),
        panel.grid.minor.y = element_blank(),
        axis.text.x        = element_text(size = 7, face = "bold"),
        legend.key.height  = grid::unit(0.28, "cm"),
        legend.key.width   = grid::unit(0.28, "cm"),
        legend.key.size    = grid::unit(0.28, "cm"),
        legend.spacing.x   = grid::unit(0.05, "cm"),
        legend.text        = element_text(size = 7, face = "bold")
      )
  }
#* PCA (square aspect + top horizontal legend)
  theme_pub_pca <- function(base_size = 12, base_family = "Arial",
                            border_linewidth = 0.8,
                            major_frac = 0.25, minor_frac = 0.25,
                            pca_border_scale = 2.0) {
    theme_pub_base(base_size, base_family,
      border_linewidth = border_linewidth * pca_border_scale,
      major_frac = major_frac, minor_frac = minor_frac,
      show_major_x = TRUE, show_major_y = TRUE,
      show_minor_x = FALSE, show_minor_y = FALSE
    ) +
      ggplot2::theme(
        aspect.ratio = 1,
        plot.margin = grid::unit(c(2, 8, 8, 8), "pt"),
        legend.position = "top",
        legend.direction = "horizontal",
        legend.box = "horizontal",
        legend.box.margin = ggplot2::margin(0, 0, 0, 0),
        legend.margin = ggplot2::margin(t = 0, r = 0, b = -3, l = 0),

        # 🔽 controls symbol size
        legend.key.width = grid::unit(0.35, "cm"),
        legend.key.height = grid::unit(0.35, "cm"),
        legend.key.size = grid::unit(0.35, "cm"),

        # 🔽 distance between symbol and text
        legend.spacing.x = grid::unit(0, "cm"),

        # 🔽 fine-tunes text alignment (0 = left, 1 = right; <0 moves text closer to key)
        legend.text.align = -10,
        legend.title = ggplot2::element_blank(),
        axis.title.y = ggplot2::element_text(
          margin = ggplot2::margin(r = 0), hjust = 0.5
        )
      )
  }
#* Helpers
#+ Scale helpers (offset axes, tight expand)
  scale_x_pub <- function(lims = NULL, breaks = ggplot2::waiver(), minor_breaks = NULL,
                          expand = c(0, 0), offset_guides = TRUE) {
    ggplot2::scale_x_continuous(
      limits = lims, breaks = breaks, minor_breaks = minor_breaks, expand = expand,
      guide = if (isTRUE(offset_guides)) ggprism::guide_prism_offset() else "none"
    )
  }
  scale_y_pub <- function(lims = NULL, breaks = ggplot2::waiver(), minor_breaks = NULL,
                          expand = c(0, 0), offset_guides = TRUE) {
    ggplot2::scale_y_continuous(
      limits = lims, breaks = breaks, minor_breaks = minor_breaks, expand = expand,
      guide = if (isTRUE(offset_guides)) ggprism::guide_prism_offset() else "none"
    )
  }
#* Palettes and Colors
  #+ Color/Fill manual scales
    scale_color_variant <- function(...) ggplot2::scale_color_manual(values = variant_colors, ...)
    scale_fill_variant <- function(...) ggplot2::scale_fill_manual(values = variant_colors, ...)
  #+ Definitions
    variant_colors <- c(
      "PTC"    = "#DF8D0A",
      "FV-PTC" = "#23744E",
      "FTC"    = "#194992"
    )
    variant_light <- c(
      "PTC"    = "#e9bb71",
      "FV-PTC" = "#80ac95",
      "FTC"    = "#7a92bc"
    )
    LVI_colors <- c(
      "-LVI" = "#9C27B0", # purple
      "+LVI" = "#FFC107" # amber
    )
    cluster_colors <- c(
      "Cluster 1" = "#94001E",
      "Cluster 2" = "#03507D"
    )
    T_stage_colors <- c(
      "Cluster 1" = "#94001E",
      "Cluster 2" = "#03507D"
    )
    cluster_light <- c(
      "Cluster 1" = "#bb697b",
      "Cluster 2" = "#6f96b1"
    )
    sex_colors <- c(
      "Male"   = "#0a9af3",
      "Female" = "#d55e70"
    )
    T_stage_bin_colors <- c(
      "T1-T2" = "#dfba37",
      "T3-T4" = "#72061c"
    )
    T_stage_colors_heatmap <- c(
      "T1" = "#f2a3b3",
      "T2" = "#f2a3b3",
      "T3" = "#72061c",
      "T4" = "#72061c"
    )
    T_stage_cluster_colors <- c(
      "T1" = "#FFF096",
      "T2" = "#FDC586",
      "T3" = "#F45158",
      "T4" = "#72061c"
    )
  #+ List of palettes
    .palettes <- list(
      variant   = variant_colors,
      variant_light = variant_light,
      LVI       = LVI_colors,
      cluster   = cluster_colors,
      cluster_light = cluster_light,
      sex       = sex_colors,
      T_bin     = T_stage_bin_colors,
      T_heatmap = T_stage_colors_heatmap,
      T_cluster = T_stage_cluster_colors
    )
  #+ Get palette helper function
    get_palette <- function(name) {
      if (!name %in% names(.palettes)) stop(sprintf("Palette '%s' not found.", name))
      .palettes[[name]]
    }
  #+ Scale helpers (one-liners you can use in plots)
    scale_color_variant <- function(...) ggplot2::scale_color_manual(values = variant_colors, ...)
    scale_fill_variant <- function(...) ggplot2::scale_fill_manual(values = variant_colors, ...)
    scale_color_LVI <- function(...) ggplot2::scale_color_manual(values = LVI_colors, ...)
    scale_fill_LVI <- function(...) ggplot2::scale_fill_manual(values = LVI_colors, ...)
    scale_color_cluster <- function(...) ggplot2::scale_color_manual(values = cluster_colors, ...)
    scale_fill_cluster <- function(...) ggplot2::scale_fill_manual(values = cluster_colors, ...)
    scale_color_sex <- function(...) ggplot2::scale_color_manual(values = sex_colors, ...)
    scale_fill_sex <- function(...) ggplot2::scale_fill_manual(values = sex_colors, ...)
    scale_color_T_bin <- function(...) ggplot2::scale_color_manual(values = T_stage_bin_colors, ...)
    scale_fill_T_bin <- function(...) ggplot2::scale_fill_manual(values = T_stage_bin_colors, ...)
    scale_color_T_cluster <- function(...) ggplot2::scale_color_manual(values = T_stage_cluster_colors, ...)
    scale_fill_T_cluster <- function(...) ggplot2::scale_fill_manual(values = T_stage_cluster_colors, ...)
    scale_fill_T_heatmap <- function(...) ggplot2::scale_fill_manual(values = T_stage_colors_heatmap, ...)
  #+ Safety check: warn if data has levels with no color =====
    check_palette_levels <- function(x, palette, palette_name = deparse(substitute(palette))) {
      lv <- if (is.factor(x)) levels(x) else sort(unique(as.character(x)))
      missing <- setdiff(lv, names(palette))
      if (length(missing)) {
        warning(sprintf(
          "Palette '%s' missing %d level(s): %s",
          palette_name, length(missing), paste(missing, collapse = ", ")
        ), call. = FALSE)
      }
      invisible(TRUE)
    }

