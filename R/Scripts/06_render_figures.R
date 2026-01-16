#* 6: Render Figures
#+ 6.1: Figure 1
fig1 <- ggdraw(xlim = c(0, 8.5), ylim = c(0, 11)) +
  # 1A
  draw_plot(p1A, x = 0.8, y = 6.683333333, width = 3.25, height = 3.1) +
  # 1B
  draw_plot(ggdraw() + draw_grob(p1B), x = 4.063333333, y = 6.62, width = 3.933333333, height = 3.6) +
  # 1C
  draw_plot(ggdraw() + draw_grob(p1C), x = 1.378333334, y = 1.913333333, width = 5.7421875, height = 5.7421875) +
  # Labels
  figure_labels(list(
    A = c(0.88, 10.04667),
    B = c(4.263333333, 10.04667),
    C = c(0.88, 6.59667),
    "Figure 1" = c(0.49, 10.43)
  ))
#+ 6.2: Figure 2
fig2 <- ggdraw(xlim = c(0, 8.5), ylim = c(0, 11)) +
  draw_plot(p2A, x = 0.75, y = 1.823333334, width = 3.75, height = 8.333333333) +
  # Manually add dagger symbol to 9-O-Acetyl-Neu5Ac
  draw_text("†", x = 1.985000001, y = 6.025, size = 4, fontface = "bold", family = "Arial") +
  # 2B
  draw_plot(p2B.1, x = 4.703333333, y = 8.456666667, width = 1.5, height = 1.5) +
  draw_plot(p2B.2, x = 6.253333333, y = 8.456666667, width = 1.5, height = 1.5) +
  draw_plot(p2B.3, x = 4.703333333, y = 7.006666667, width = 1.5, height = 1.5) +
  draw_plot(p2B.4, x = 6.253333333, y = 7.006666667, width = 1.5, height = 1.5) +
  draw_plot(p2B.5, x = 4.703333333, y = 5.556666667, width = 1.5, height = 1.5) +
  draw_plot(p2B.6, x = 6.253333333, y = 5.556666667, width = 1.5, height = 1.5) +
  # 2C
  draw_plot(p2C.1, x = 4.703333333, y = 3.77, width = 1.5, height = 1.5) +
  draw_plot(p2C.2, x = 6.253333333, y = 3.77, width = 1.5, height = 1.5) +
  draw_plot(p2C.3, x = 4.703333333, y = 2.32, width = 1.5, height = 1.5) +
  draw_plot(p2C.4, x = 6.253333333, y = 2.32, width = 1.5, height = 1.5) +
  # Labels
  figure_labels(list(
    A = c(0.8266666667, 10.04667),
    B = c(4.753333333, 10.04667),
    C = c(4.746666666, 5.36667),
    "Figure 2" = c(0.49, 10.43)
  ))
#+ 6.3: Figure 3
fig3 <- ggdraw(xlim = c(0, 8.5), ylim = c(0, 11)) +
  # 3A/C
  draw_plot(ggdraw() + draw_grob(p3AC), x = 0.6033333334, y = 1.25, width = 6.885, height = 8.91)+
  # 3B
  draw_plot(p3B.R1.C1 , x = 5.04, y = 8.403333334, width = 1.5, height = 1.5) +
  draw_plot(p3B.R1.C2, x = 6.479999999, y = 8.403333334, width = 1.5, height = 1.5) +
  draw_plot(p3B.R2.C1 , x = 5.04, y = 6.97, width = 1.5, height = 1.5) +
  draw_plot(p3B.R2.C2, x = 6.479999999, y = 6.97, width = 1.5, height = 1.5) +
  # 3D
  draw_plot(p3D.R1.C1 , x = 5.04, y = 5.124993367, width = 1.5, height = 1.5) +
  draw_plot(p3D.R1.C2, x = 6.479999999, y = 5.124993367, width = 1.5, height = 1.5) +
  draw_plot(p3D.R2.C1 , x = 5.04, y = 3.691660033, width = 1.5, height = 1.5) +
  draw_plot(p3D.R2.C2, x = 6.479999999, y = 3.691660033, width = 1.5, height = 1.5) +
  draw_plot(p3D.R3.C1 , x = 5.04, y = 2.191660033, width = 1.5, height = 1.5) +
  draw_plot(p3D.R3.C2, x = 6.479999999, y = 2.191660033, width = 1.5, height = 1.5) +
  # Labels
  figure_labels(list(
    A = c(0.9900003333, 10.04667),
    B = c(5.13, 10.04667),
    C = c(0.9900003, 6.94),
    D = c(5.13, 6.94),
    "Figure 3" = c(0.49, 10.43)
  ))
#+ 6.4: Figure 4
fig4 <- ggdraw(xlim = c(0, 8.5), ylim = c(0, 11)) +
  draw_plot(ggdraw() + draw_grob(p4), x = 0, y = -0.295, width = 8.5, height = 11) +
  # Labels
  figure_labels(list(
    "Figure 4" = c(0.49, 10.43)
  ))
#+ 6.5: Print all figures
#- 6.5.1: Print PNGs
print_to_png(fig1, "Final/PNG/fig1.png")
print_to_png(fig2, "Final/PNG/fig2.png")
print_to_png(fig3, "Final/PNG/fig3.png")
print_to_png(fig4, "Final/PNG/fig4.png")
#- 6.5.2: Print TIFFs
print_to_tiff(fig1, "Final/TIFF/fig1.tiff")
print_to_tiff(fig2, "Final/TIFF/fig2.tiff")
print_to_tiff(fig3, "Final/TIFF/fig3.tiff")
print_to_tiff(fig4, "Final/TIFF/fig4.tiff")
#- 6.5.3: Print PDFs
{
  # Figure 1
  tryCatch({
    pdf("Outputs/Figures/Final/PDF/Figure 1.pdf", width = 8.5, height = 11)
    img1 <- readPNG("Outputs/Figures/Final/PNG/fig1.png", native = TRUE)
    grid.newpage()
    grid.raster(img1, width = unit(8.5, "inches"), height = unit(11, "inches"), interpolate = TRUE)
    dev.off()
  }, error = function(e) {
    dev.off()
    cat("Warning: Figure 1 PDF failed:", e$message, "\n")
  })
  # Figure 2
  tryCatch({
    pdf("Outputs/Figures/Final/PDF/Figure 2.pdf", width = 8.5, height = 11)
    img2 <- readPNG("Outputs/Figures/Final/PNG/fig2.png", native = TRUE)
    grid.newpage()
    grid.raster(img2, width = unit(8.5, "inches"), height = unit(11, "inches"), interpolate = TRUE)
    dev.off()
  }, error = function(e) {
    dev.off()
    cat("Warning: Figure 2 PDF failed:", e$message, "\n")
  })
  
  # Figure 3
  tryCatch({
    pdf("Outputs/Figures/Final/PDF/Figure 3.pdf", width = 8.5, height = 11)
    img3 <- readPNG("Outputs/Figures/Final/PNG/fig3.png", native = TRUE)
    grid.newpage()
    grid.raster(img3, width = unit(8.5, "inches"), height = unit(11, "inches"), interpolate = TRUE)
    dev.off()
  }, error = function(e) {
    dev.off()
    cat("Warning: Figure 3 PDF failed:", e$message, "\n")
  })
  # Figure 4
  tryCatch({
    pdf("Outputs/Figures/Final/PDF/Figure 4.pdf", width = 8.5, height = 11)
    img4 <- readPNG("Outputs/Figures/Final/PNG/fig4.png", native = TRUE)
    grid.newpage()
    grid.raster(img4, width = unit(8.5, "inches"), height = unit(11, "inches"), interpolate = TRUE)
    dev.off()
  }, error = function(e) {
    dev.off()
    cat("Warning: Figure 4 PDF failed:", e$message, "\n")
  })
  cat("✓ PDF generation complete\n")
}
#+ 6.6: Create compiled PDF
pdf_combine(
  input = c(
    "Outputs/Figures/Final/PDF/Figure 1.pdf",
    "Outputs/Figures/Final/PDF/Figure 2.pdf",
    "Outputs/Figures/Final/PDF/Figure 3.pdf",
    "Outputs/Figures/Final/PDF/Figure 4.pdf"
  ),
  output = "Outputs/Figures/Final/Figs1-4.pdf"
)
