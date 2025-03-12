# testing, https://ggplot2-book.org/ext-springs#testing-the-stat

library(ggplot2)
library(tibble)
library(dplyr)

set.seed(421)
df <- tibble(
  x = runif(5, max = 10),
  y = runif(5, max = 10),
  xend = runif(5, max = 10),
  yend = runif(5, max = 10),
  class = sample(letters[1:2], 5, replace = TRUE)
)

ggplot(df) +
  geom_spring(aes(x = x, y = y,
                  xend = xend, yend = yend,
                  color = class),
              linewidth = 2) +
  coord_equal()

# faceting

ggplot(df) +
  geom_spring(
    aes(x, y, xend = xend, yend = yend, color = class),
    linewidth = 1.5
  ) +
  coord_equal() +
  facet_wrap(~ class)

# tension and diameter as aesthetics


df <- tibble(
  x = runif(5, max = 10),
  y = runif(5, max = 10),
  xend = runif(5, max = 10),
  yend = runif(5, max = 10),
  class = sample(letters[1:2], 5, replace = TRUE),
  tension = runif(5),
  diameter = runif(5, 0.5, 1.5)
)

ggplot(df, aes(x, y, xend = xend, yend = yend)) +
  geom_spring(aes(tension = tension,
                  diameter = diameter,
                  color = class),
              linewidth = 1.2)

