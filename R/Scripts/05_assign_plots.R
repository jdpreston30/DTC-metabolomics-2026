#* 5: Assign and Render Plots
#+ 5.1: Figure 1
p1A <- volcano
p1B <- grid::rasterGrob(as.raster(magick::image_read("Outputs/Figures/Raw/mfn_enrich.png")))
p1C <- grid::rasterGrob(as.raster(magick::image_read("Outputs/Figures/Raw/mfn_network_laptop.png")))
#+ 5.2: Figure 2
p2A <- div_bars
p2B.1 <- stage_feature_plots$GMP
p2B.2 <- stage_feature_plots$AMP
p2B.3 <- stage_feature_plots$Oleate
p2B.4 <- stage_feature_plots$`γ-Linolenate`
p2B.5 <- stage_feature_plots$SAH
p2B.6 <- stage_feature_plots$Kynurenine
p2C.1 <- stage_feature_plots$`2,3-Dihydroxybenzoate`
p2C.2 <- stage_feature_plots$`α-Ketoisocaproate`
p2C.3 <- stage_feature_plots$`Acetyl phosphate`
p2C.4 <- stage_feature_plots$Adrenaline
#+ 5.3: Figure 3
p3A_legend <- plot_corr_legend()
p3A <- grid::rasterGrob(as.raster(
  magick::image_read("Outputs/Figures/Raw/p3A.png") %>%
  magick::image_crop("5412x4601+0+1065")))
#+ 5.5: Supplementary 1
kegg_enrichment_plot
