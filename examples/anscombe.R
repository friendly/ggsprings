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

library(ggplot2)
library(dplyr)
data(anscombe)

theme_set(theme_minimal(base_size = 16))

simple_plot <-
  ggplot(anscombe,
         aes(x = x1, y = y1)) +
  geom_point(size = 4, color = "red") +
  geom_smooth(method = lm, formula = y~x,
              se = FALSE,
              linewidth = 1.4)
simple_plot

model <- lm(y1 ~ x1, data = anscombe)

anscombe <- anscombe |>
  mutate(fitted = predict(model))

# doesn't work this way because the original anscombe w/o fitted was captured in simple plot
simple_plot +
  geom_spring(aes(xend = x1,
                  yend = fitted),
              color = "blue",
              linewidth = 1)

# supply new data arg to geom_spring
simple_plot +
  geom_spring(data=anscombe,
              aes(xend = x1, yend = fitted),
              color = "blue",
              linewidth = 1)

# Using broom::augment
broom::augment(model)

ans_mod <- lm(y1 ~ x1, data = anscombe)
ans_aug <- broom::augment(ans_mod)

ggplot(ans_aug) +
  aes(x = x1, y = y1) +
  geom_point(size = 4, color = "red") +
  geom_smooth(method = lm, formula = y~x,
              se = FALSE,
              linewidth = 1.4) +
  geom_spring(
            aes(xend = x1, yend = .fitted),
            color = "blue",
            linewidth = 1)


# x2, y2


ggplot(anscombe,
       aes(x = x2, y = y2)) +
  geom_point(size = 4, color = "red") +
  geom_smooth(method = lm, formula = y~x,
              se = FALSE,
              linewidth = 1.4)


  # ggplot(data = anscombe, aes(x = x1, y = y1)) +
  # geom_point(size = 5, color = "red") +
  # stat_smooth(method = "lm", se = FALSE) +
  # geom_spring(aes(xend = x1,
  #                 yend = fitted),
  #             color = "blue",
  #             linewidth = 1) +
  # theme_bw(base_size = 14)


