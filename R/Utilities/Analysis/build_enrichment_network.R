# EnrichNet-style plot with:
# - node color = p-value (blues, dark=small p, light=large p)
# - solid filled circles with black outlines
# - all text in black Arial
# - optional p-value limits for the colorbar (default 0.01–0.05 as you showed)
# Drop-in replacement for the plotting tail of build_enrichment_network(),
# or use this full function if you want it turnkey.

build_enrichment_network <- function(
    enrich_df, # cols: pathway_ID, pathway_name, enrichment_factor, neg_log_p OR p_value
    edge_thresh = 0.10,
    prefer_hsa = TRUE,
    term2compound_override = NULL,
    save_path = NULL,
    width = 8, height = 6, dpi = 300, units = "in", bg = "white",
    seed = 123, layout = "fr",
    p_limits = c(0.01, 0.05) # <- your preferred color scale limits
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

  # ---- PLOT: blue gradient, solid filled circles, black outlines, Arial text ----
  p <- ggraph::ggraph(g, layout = layout) +
    ggraph::geom_edge_link(aes(width = weight), alpha = 0.5, colour = "grey50") +

    # Nodes: size = enrichment factor, fill = p-value color
    ggraph::geom_node_point(aes(size = size_val, fill = p_val),
      shape = 21, stroke = 0.6, colour = "black"
    ) +

    # Labels: bold, below the dots
    ggraph::geom_node_text(aes(label = Name),
      size = 4, colour = "black",
      family = "Arial", fontface = "bold",
      vjust = 2.5
    ) +

    # Enrichment factor (size legend)
    scale_size_continuous(
      range = c(3,30),
      limits = c(0, 5),
      breaks = c(4, 3, 2, 1), # show big → small
      name = "\nEnrichment factor",
      guide = guide_legend(
        reverse = TRUE, # reverse so 4 is at top
        override.aes = list(fill = "black", colour = "black")
      )
    ) +

    # p-value (fill legend)
    scale_fill_gradient(
      low = "#0a2256", high = "#c3dbe9",
      limits = c(0.01, 0.05), oob = scales::squish,
      name = "p-value\n",
      guide = guide_colorbar(
        reverse = TRUE, # 0.01 at top, 0.05 bottom
        barheight = grid::unit(5, "cm"),
        barwidth = grid::unit(0.9, "cm")
      )
    ) +

    ggraph::scale_edge_width(range = c(0.2, 2), guide = "none") +
    theme_void(base_family = "Arial") +
    theme(
      legend.position = "right",
      legend.margin   = margin(l = 50, r = 10, t = 0, b = 0),
      legend.text     = element_text(size = 12, face = "bold", family = "Arial"), # match node labels
      legend.title    = element_text(size = 13, face = "bold", family = "Arial"),
      text            = element_text(color = "black", family = "Arial")
      # plot.margin     = margin(5, 5, 5, 5)
    ) +
  ggplot2::coord_equal(clip = "off") +
  ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = 0.15)) + # left/right padding
  ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = 0.25))
  #- Saving
  if (!is.null(save_path)) {
    ggplot2::ggsave(
      filename = save_path, plot = p,
      width = width, height = height, dpi = dpi, units = units, bg = bg
    )
  }

  list(plot = p, graph = g, nodes = nodes, edges = edges, term2compound = t2c)
}
# keep everything visible beyond the panel; add generous padding
