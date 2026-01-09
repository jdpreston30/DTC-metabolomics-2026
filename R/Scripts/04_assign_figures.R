#* 4: Assign and Render Plots
#+ 4.1: Assign Plots
A <- volcano
B <- grid::rasterGrob(as.raster(magick::image_read("Outputs/Figures/Raw/mfn_enrich.png")))
C <- grid::rasterGrob(as.raster(magick::image_read("Outputs/Figures/Raw/mfn_network.png")))
D.1 <- stage_feature_plots$GMP + 
  theme(
    axis.text.x = element_text(color = "transparent"),
    axis.ticks.x = element_line(color = "transparent")
  )
D.2 <- stage_feature_plots$AMP + 
  theme(
    axis.text.x = element_text(color = "transparent"),
    axis.ticks.x = element_line(color = "transparent"),
    axis.title.y = element_text(color = "transparent")
  )
D.3 <- stage_feature_plots$`γ-Linolenic Acid`
D.4 <- stage_feature_plots$SAH + 
  theme(
    axis.title.y = element_text(color = "transparent")
  )
