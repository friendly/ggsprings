# for an icon
library(ggplot2)
library(tibble)

df <- tibble(
  x = 1:4,
  y = x/2 + c(3, -1.5, 3, -2.5)
)

mod <- lm(y ~ x, data = df)

df_aug <- broom::augment(mod)

ggplot(df_aug, aes(x=x, y=y)) +
  geom_point(size=3) +
  geom_smooth(method = lm, se = FALSE)


