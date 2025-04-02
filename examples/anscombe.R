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

# ----------

anscombe |>
  ggplot() +
  aes(x = x1, y = y1) +
  geom_point(size = 5, color = "red") +
  geom_smooth(method = lm) +
  theme_minimal(base_size = 15)


library(dplyr)
model <- lm(y1 ~ x1, data = anscombe)

anscombe <- anscombe |>
  mutate(fitted = predict(model))


ggplot(data = anscombe, aes(x = x1, y = y1)) +
  geom_point(size = 5, color = "red") +
  stat_smooth(method = "lm", se = FALSE) +
  geom_spring(aes(xend = x1,
                  yend = fitted),
              color = "blue",
              linewidth = 1) +
  theme_bw(base_size = 14)


  # geom_spring(aes(x = gdpPercap,
  #                 xend = gdpPercap,
  #                 y = lifeExp,
  #                 yend = yhat,
  #                 diameter = diameter,
  #                 tension = tension), color = "darkgray") +



