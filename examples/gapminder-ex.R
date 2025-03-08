library(ggplot2)
library(dplyr)


source(here::here("examples/springs.R"))

data(gapminder, package = "gapminder")
gm2007 <- gapminder |>
  filter(year == 2007, continent == "Asia")

regression_model <- lm(lifeExp ~ log10(gdpPercap), gm2007)


# get fitted values, and set diameters and tension for springs
gm2007 <- gm2007 |>
  mutate(yhat = predict(regression_model),
          diameter = .02,
          tension = 5 + (lifeExp - yhat)^2)

# means of X and Y
xbar <- 10^mean(log10(gm2007$gdpPercap))
ybar <- mean(gm2007$lifeExp)

# plot gdpPercap vs. lifeExp
simple_plot <- ggplot(gm2007, aes(x = gdpPercap, y = lifeExp)) +
  geom_point(size = 2) +
  geom_point(aes(x = xbar, y = ybar), color = "blue", size = 3) +
  scale_x_log10() +
  ylim(c(40, 85))

# add springs
spring_plot <- simple_plot +
  geom_spring(aes(x = gdpPercap,
                  xend = gdpPercap,
                  y = lifeExp,
                  yend = yhat,
                  diameter = diameter,
                  tension = tension), color = "darkgray") +
  stat_smooth(method = "lm", se = FALSE) +
  geom_point(size = 2)

spring_plot


