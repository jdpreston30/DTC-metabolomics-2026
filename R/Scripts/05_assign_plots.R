#* 5: Assign and Render Plots
#+ 5.1: Figure 1
p1A <- volcano
p1B <- grid::rasterGrob(as.raster(magick::image_read("Outputs/Figures/Raw/p1B.png")))
p1C <- grid::rasterGrob(as.raster(magick::image_read("Outputs/Figures/Raw/p1C.png")))
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
#- 5.3.1: Read in raw for A/C
p3A_legend <- plot_corr_legend()
p3A <- grid::rasterGrob(as.raster(
  magick::image_read("Outputs/Figures/Raw/p3A.png") %>%
  magick::image_crop("5412x4601+0+1065")))
p3C_legend <- plot_corr_mummi_legend()
p3C <- grid::rasterGrob(as.raster(
  magick::image_read("Outputs/Figures/Raw/p3C.png") %>%
  magick::image_crop("5800x7000+900+0")))
#- 5.3.2: Knit together A/C
fig_3AC <- ggdraw(xlim = c(0, 8.5), ylim = c(0, 11)) +
  draw_plot(ggdraw() + draw_grob(p3A), x = 0.473333333, y = 7.470000001, width = 3.996, height = 3.6) +
  draw_plot(ggdraw() + draw_grob(p3C), x = -0.02, y = 1.836666667, width = 4.511700458, height = 5.445129867) +
  draw_plot(p3A_legend, x = 4.538333, y = 6.948333334, width = 0.9, height = 4.185) +
  draw_plot(p3C_legend, x = 4.501667, y = 1.975000000, width = 0.9, height = 4.185)
#- 5.3.3: Save knit A/C
print_to_png(fig_3AC, "Raw/fig_3AC.png", dpi = 1200, background = "transparent")
#- 5.3.4: Read back in and crop
p3AC <- grid::rasterGrob(as.raster(magick::image_read("Outputs/Figures/Raw/fig_3AC.png")))
#- 5.3.5: Assign 3B plots
p3B.R1.C1 <- AP_GMP
p3B.R1.C2 <- Adr_3KS
p3B.R2.C1 <- Adr_Kyn
p3B.R2.C2 <- Adr_SAH
#- 5.3.6: Assign 3D plots
p3D.R1.C1 <- GMP_R5P
p3D.R1.C2 <- KS3_ODHAP
p3D.R2.C1 <- MNA1_SAH
p3D.R2.C2 <- PAPS_SAH
p3D.R3.C1 <- AcGlu_Cit
p3D.R3.C2 <- Kyn_Ser