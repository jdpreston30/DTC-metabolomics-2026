create_4up_preview <- function(plot, 
                                plot_width = 3, 
                                plot_height = 2.5,
                                plot_margin = 4,
                                ncol = 2, 
                                nrow = 2,
                                page_margins = c(1, 1, 1, 1),
                                title = "Figure 3",
                                tag_levels = "A",
                                collect_guides = FALSE) {
  
  # Create individual plots with controlled dimensions
  p_small <- plot + 
    ggplot2::theme(
      plot.margin = grid::unit(rep(plot_margin, 4), "pt"),
      aspect.ratio = plot_height/plot_width
    )

  # Create layout with precise control
  page_4up <- patchwork::wrap_plots(
    p_small, p_small,
    p_small, p_small,
    ncol = ncol, 
    nrow = nrow,
    widths = rep(plot_width, ncol),
    heights = rep(plot_height, nrow),
    guides = if(collect_guides) "collect" else "auto"
  ) +
    patchwork::plot_annotation(
      title = title,
      tag_levels = tag_levels,
      theme = ggplot2::theme(
        plot.title.position = "plot",
        plot.title = ggplot2::element_text(hjust = 0, face = "bold", family = "Arial", size = 15),
        plot.margin = grid::unit(page_margins, "in")
      )
    )

  smart_preview(page_4up)
  return(page_4up)
}
