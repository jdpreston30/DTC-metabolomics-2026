  s <- -4.14
fig3 <- ggdraw(xlim = c(0, 8.5), ylim = c(0, 11)) +
  # 3A
  draw_plot(ggdraw() + draw_grob(p3A_combined), x = 0.6033333334, y = 1.25, width = 6.885, height = 8.91) +
  # 3B
  draw_plot(SAM_SAH_cor , x = 5.04, y = 8.403333334, width = 1.5, height = 1.5) +
  # draw_plot(p3B.R1.C2, x = 6.479999999, y = 8.403333334, width = 1.5, height = 1.5) +
  # draw_plot(p3B.R2.C1 , x = 5.04, y = 6.97, width = 1.5, height = 1.5) +
  # draw_plot(p3B.R2.C2, x = 6.479999999, y = 6.97, width = 1.5, height = 1.5) +
######## TEMPORARY TEST ########
  # 3B
  draw_plot(p3B.R1.C1 , x = 5.04+s, y = 5.124993367, width = 1.5, height = 1.5) +
  draw_plot(p3B.R1.C2, x = 6.479999999+s, y = 5.124993367, width = 1.5, height = 1.5) +
  draw_plot(p3B.R2.C1 , x = 5.04+s, y = 3.691660033, width = 1.5, height = 1.5) +
  draw_plot(p3B.R2.C2, x = 6.479999999+s, y = 3.691660033, width = 1.5, height = 1.5) +
  # 3C
  draw_plot(p3C.R1.C1 , x = 5.04, y = 5.124993367, width = 1.5, height = 1.5) +
  draw_plot(p3C.R1.C2, x = 6.479999999, y = 5.124993367, width = 1.5, height = 1.5) +
  draw_plot(p3C.R2.C1 , x = 5.04, y = 3.691660033, width = 1.5, height = 1.5) +
  draw_plot(p3C.R2.C2, x = 6.479999999, y = 3.691660033, width = 1.5, height = 1.5) +
  draw_plot(p3C.R3.C1 , x = 5.04, y = 2.191660033, width = 1.5, height = 1.5) +
  draw_plot(p3C.R3.C2, x = 6.479999999, y = 2.191660033, width = 1.5, height = 1.5) +
  # Labels
  figure_labels(list(
    A = c(0.9900003333, 10.04667),
    B = c(5.13+s, 6.94),
    C = c(5.13, 6.94),
    "Figure 3" = c(0.49, 10.43)
  ))
  print_to_png(fig3, "fig3_preview.png", dpi = 300)