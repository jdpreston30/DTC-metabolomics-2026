#* 4: Assign and Render Plots
#+ 4.1: Figure 1
p1A <- volcano
p1B <- grid::rasterGrob(as.raster(magick::image_read("Outputs/Figures/Raw/mfn_enrich.png")))
p1C <- grid::rasterGrob(as.raster(magick::image_read("Outputs/Figures/Raw/mfn_network_laptop.png")))
#+ 4.1: Figure 2
# p2A.1 <- stage_feature_plots$GMP + 
#   theme(
#     axis.text.x = element_text(color = "transparent"),
#     axis.ticks.x = element_line(color = "transparent")
#   )
# p2A.2 <- stage_feature_plots$AMP + 
#   theme(
#     axis.text.x = element_text(color = "transparent"),
#     axis.ticks.x = element_line(color = "transparent"),
#     axis.title.y = element_text(color = "transparent")
#   )
# p2A.3 <- stage_feature_plots$`γ-Linolenic Acid`
# p2A.4 <- stage_feature_plots$SAH + 
  theme(
    axis.title.y = element_text(color = "transparent")
  )
