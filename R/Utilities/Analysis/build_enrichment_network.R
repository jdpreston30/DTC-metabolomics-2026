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
