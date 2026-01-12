#* 5: Render Figures
#+ 5.1: Figure 1
fig1 <- ggdraw(xlim = c(0, 8.5), ylim = c(0, 11)) +
  draw_plot(p1A, x = 1.0, y = 6.60, width = 3.25, height = 3.1) +
  draw_plot(ggdraw() + draw_grob(p1B), x = 4.2, y = 6.62, width = 5.9 / 1.5, height = 5.4 / 1.5) +
  draw_plot(ggdraw() + draw_grob(p1C), x = 1, y = 3.7, width = 10.5 / 3.2, height = 10.5 / 3.2) +
  figure_labels(list(
    A = c(1.08, 9.88),
    B = c(4.40, 9.88),
    C = c(1.08, 6.43)
  ))
#+ 5.2: Figure 2
fig2 <- ggdraw(xlim = c(0, 8.5), ylim = c(0, 11)) +
  draw_plot(D.1, x = 4.75, y = 5, width = 1.5, height = 1.5) +
  draw_plot(D.2, x = 6.2, y = 5, width = 1.5, height = 1.5) +
  draw_plot(D.3, x = 4.75, y = 3.65, width = 1.5, height = 1.5) +
  draw_plot(D.4, x = 6.2, y = 3.65, width = 1.5, height = 1.5) +
  figure_labels(list(
    D = c(4.40, 6.43)
  ))
#+ 5.3: Print all figures
print_to_tiff(fig1, "Final/fig1.tiff")
print_to_tiff(fig2, "Final/fig2.tiff")