# illustrate mean by springs

library(tibble)

set.seed(1234)
N <- 8
df <- tibble(
  x = runif(N, 1, 10),
  y = seq(1, N)
)

means <- colMeans(df)
xbar <- means[1] |> print()

df <- df |>
  mutate(tension = (x - xbar)^2,
         diameter = 0.4)

ggplot(df, aes(x=x, y=y)) +
  geom_point(size = 5, color = "red") +
  geom_segment(x = xbar, xend = xbar,
               y = 1/2,  yend = N + 1/2,
               linewidth = 3) +
  geom_spring(aes(x = x, xend = xbar,
                  y = y, yend = y,
                  tension = tension / 4,
                  diameter = diameter),
              color = "blue",
              linewidth = 1.2) +
  labs(x = "Value (x)",
       y = "Observation number") +
  ylim(0, N+2) +
  scale_y_continuous(breaks = 1:N) +
  annotate("text", x = xbar, y = N + 1,
           label = "Moveable\nrod", size = 5,
           , lineheight = 3/4) +
  theme_minimal(base_size = 15)


# centroid by springs

set.seed(1234)
N <- 10
df <- tibble(
  x = runif(N, 0, 10),
  y = runif(N, 0, 10)
)

means <- colMeans(df)
xbar <- means[1]; ybar <- means[2]

df <- df |>
  mutate(tension = (x - xbar)^2 + (y - ybar)^2,
         diameter = 0.4)


ggplot(df, aes(x=x, y=y)) +
  geom_point(size = 5, color = "red") +
  geom_spring(aes(x = x, xend = xbar,
                  y = y, yend = ybar,
                  tension = tension / 5,
                  diameter = diameter),
              color = "blue",
              linewidth = 1.2) +
  geom_point(x = xbar, y = ybar,
             size = 7,
             shape = 15,
             color = "black") +
  # scale_x_continuous(breaks = 1:10) +
  # scale_y_continuous(breaks = 1:10) +
#  xlim(1, 10) + ylim(1, 10) +
  theme_minimal(base_size = 15)


