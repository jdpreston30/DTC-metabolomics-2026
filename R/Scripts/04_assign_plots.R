#* 4: Assign and Render Plots
#+ 4.1: Figure 1
p1A <- mad_500$heatmap_plot$gtable
p1B <- volcano
p1C <- grid::rasterGrob(as.raster(magick::image_read("Outputs/Figures/Raw/p1C.png")))
p1D <- grid::rasterGrob(as.raster(magick::image_read("Outputs/Figures/Raw/p1D.png")))
#+ 4.2: Figure 2
p2A <- div_bars
p2B.1 <- stage_feature_plots$GMP
p2B.2 <- stage_feature_plots$AMP
p2B.3 <- stage_feature_plots$Oleate
p2B.4 <- stage_feature_plots$`γ-Linolenate`
p2B.5 <- stage_feature_plots$SAH
p2B.6 <- stage_feature_plots$`Kynurenine*`
p2C.1 <- stage_feature_plots$`2,3-Dihydroxybenzoate`
p2C.2 <- stage_feature_plots$`α-Ketoisocaproate`
p2C.3 <- stage_feature_plots$`Acetyl phosphate`
p2C.4 <- stage_feature_plots$Adrenaline
#+ 4.3: Figure 3
#- 4.3.1: Read in raw for A/C
p3A_legend <- plot_corr_legend()
p3A <- grid::rasterGrob(as.raster(
  magick::image_read("Outputs/Figures/Raw/p3A.png") %>%
  magick::image_crop("5412x4601+0+1065")))
#- 4.3.2: Knit together A plot and legend
fig_3A <- ggdraw(xlim = c(0, 8.5), ylim = c(0, 11)) +
  draw_plot(ggdraw() + draw_grob(p3A), x = 0.473333333, y = 7.470000001, width = 3.996, height = 3.6) +
  draw_plot(p3A_legend, x = 4.538333, y = 6.948333334, width = 0.9, height = 4.185)
#- 4.3.3: Save knit A/C
print_to_png(fig_3A, "Raw/fig_3A.png", dpi = 1200, background = "transparent")
#- 4.3.4: Read back in and crop
p3A_combined <- grid::rasterGrob(as.raster(magick::image_read("Outputs/Figures/Raw/fig_3A.png")))
#- 4.3.5: Assign 3B plots
p3B.R1.C1 <- AP_GMP
p3B.R1.C2 <- Adr_3KS
p3B.R2.C1 <- Adr_Kyn
p3B.R2.C2 <- Adr_SAH
#- 4.3.6: Assign 3C plots
p3C.R1.C1 <- GMP_R5P
p3C.R1.C2 <- KS3_ODHAP
p3C.R1.C3 <- MNA1_SAH
p3C.R2.C1 <- PAPS_SAH
p3C.R2.C2 <- AcGlu_Cit
p3C.R2.C3 <- Kyn_Ser
#- 4.3.7: Assign 3D plots
p3D.R1.C1 <- ratio_MNA1
p3D.R1.C2 <- ratio_gLin
p3D.R2.C1 <- ratio_GMP
p3D.R2.C2 <- ratio_AMP
#+ 4.4: Figure 4
p4 <- grid::rasterGrob(as.raster(magick::image_read("Outputs/Figures/Raw/p4.png")))
