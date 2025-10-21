#* 4: Assign and Render Plots
#+ 4.1: Assign Plots
volcano <- volcano
enrich <- grid::rasterGrob(as.raster(magick::image_read("Outputs/Figures/Raw/mfn_enrich.png")))
network <- grid::rasterGrob(as.raster(magick::image_read("Outputs/Figures/Raw/mfn_network.png")))
#+ 4.2: Render Figures
print_to_png(ggdraw(xlim = c(0, 8.5), ylim = c(0, 11)) +
  #- 2A
  draw_plot(volcano, x = .93, y = 6.95, width = 3.25, height = 3.1) +
  #- 2B
  draw_plot(ggdraw() + draw_grob(enrich), x = 4.15, y = 7.27, width = 5.9 / 1.7, height = 5.4 / 1.7) +
  #- 2C
  #- 2D
  draw_plot(ggdraw() + draw_grob(network), x = 0.6, y = 3.9, width = 10.5 / 2.9, height = 10.5 / 2.9) + grdgd()+
  #- Labels
  figure_labels(list(
    A = c(0.83, 10.13),
    B = c(4.15, 10.13),
    C = c(0.83, 6.68),
    D = c(4.15, 6.68)
  )), "fig1.png", width = 8.5, height = 11, dpi = 200)
