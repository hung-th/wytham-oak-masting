library(ggplot2)

proj_theme <- theme_classic(base_size = 11) +
  theme(
    axis.title        = element_text(size = 12, colour = "black"),
    axis.text         = element_text(size = 10,  colour = "black"),
    axis.line         = element_line(linewidth = 0.6, colour = "black"),
    axis.ticks        = element_line(linewidth = 0.6, colour = "black"),
    legend.title      = element_text(size = 10,  colour = "black"),
    legend.text       = element_text(size = 9,   colour = "black"),
    strip.text        = element_text(size = 11,  colour = "black"),
    strip.background  = element_blank(),
    panel.grid        = element_blank()
  )

theme_set(proj_theme)

# NEJM palette (use as needed)
nejm_4 <- c("2020" = "#0072B5", "2021" = "#BC3C29",
             "2022" = "#E18727", "2023" = "#20854E")
nejm_3 <- c("2021" = "#BC3C29", "2022" = "#0072B5", "2023" = "#20854E")
