# Anscombe example, from
# https://github.com/EvaMaeRey/ggsmoothfit?tab=readme-ov-file#via-friendly-ggplot2-extenders-ggsprings-discussion-and-springs-extension-case-study

anscombe |>
  ggplot() +
  aes(x = x1, y = y1) +
  geom_point() +
  geom_smooth(method = lm) +
  stat_smooth_fit(geom = GeomSmoothSpring, method = lm)

last_plot() +
  aes(x = x2, y = y2)

# etc
