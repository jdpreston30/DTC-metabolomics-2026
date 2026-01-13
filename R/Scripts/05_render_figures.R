#* 5: Render Figures
#+ 5.1: Figure 1
fig1 <- ggdraw(xlim = c(0, 8.5), ylim = c(0, 11)) +
  draw_plot(p1A, x = 0.8, y = 6.683333333, width = 3.25, height = 3.1) +
  draw_plot(ggdraw() + draw_grob(p1B), x = 4.063333333, y = 6.62, width = 3.933333333, height = 3.6) +
  draw_plot(ggdraw() + draw_grob(p1C), x = 1.378333334, y = 1.913333333, width = 5.7421875, height = 5.7421875) +
  figure_labels(list(
    A = c(0.88, 9.88),
    B = c(4.263333333, 9.88),
    C = c(0.88, 6.43),
    "Figure 1" = c(0.49, 10.43)
  ))
#+ 5.2: Figure 2
fig2 <- ggdraw(xlim = c(0, 8.5), ylim = c(0, 11)) +
  draw_plot(p2A, x = 0.75, y = 1.823333334, width = 3.75, height = 8.333333333) +
  draw_plot(p2B.1, x = 4.703333333, y = 8.456666667, width = 1.5, height = 1.5) +
  draw_plot(p2B.2, x = 6.253333333, y = 8.456666667, width = 1.5, height = 1.5) +
  draw_plot(p2B.3, x = 4.703333333, y = 7.006666667, width = 1.5, height = 1.5) +
  draw_plot(p2B.4, x = 6.253333333, y = 7.006666667, width = 1.5, height = 1.5) +
  draw_plot(p2B.5, x = 4.703333333, y = 5.556666667, width = 1.5, height = 1.5) +
  draw_plot(p2B.6, x = 6.253333333, y = 5.556666667, width = 1.5, height = 1.5) +
  draw_plot(p2C.1, x = 4.703333333, y = 3.77, width = 1.5, height = 1.5) +
  draw_plot(p2C.2, x = 6.253333333, y = 3.77, width = 1.5, height = 1.5) +
  draw_plot(p2C.3, x = 4.703333333, y = 2.32, width = 1.5, height = 1.5) +
  draw_plot(p2C.4, x = 6.253333333, y = 2.32, width = 1.5, height = 1.5) +
  figure_labels(list(
    A = c(0.8266666667, 9.88),
    B = c(4.753333333, 9.88),
    C = c(4.746666666, 5.2),
    "Figure 2" = c(0.49, 10.43)
  ))
#+ 5.3: Figure 3

#+ 5.4: Print all figures
print_to_tiff(fig1, "Final/fig1.tiff")
print_to_tiff(fig2, "Final/fig2.tiff")
