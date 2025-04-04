## ----nomessages, echo = FALSE-------------------------------------------------
knitr::opts_chunk$set(
  warning = FALSE,
  message = FALSE,
  fig.height = 5,
  fig.width = 5
)
options(digits=4)
par(mar=c(3,3,1,1)+.1)
options(tinytable_html_mathjax = TRUE)

## ----setup, include=FALSE-----------------------------------------------------
library(ggsprings)
library(tibble)
library(tinytable)
library(dplyr)
library(matlib)
# set basic theme
theme_set(theme_minimal(base_size = 16))

## -----------------------------------------------------------------------------
knitr::include_graphics(here::here("man/figures/horizontal-spring2.png"))

## -----------------------------------------------------------------------------
knitr::include_graphics(here::here("man/figures/Gaines-fig2.png"))

## -----------------------------------------------------------------------------
x <- 0:6
y <- c(9, 6, 8, 4, 4, 1, 3)
lm(y~x) |> coef()

## -----------------------------------------------------------------------------
knitr::include_graphics(here::here("man/figures/Hookes-law-springs.png"))

## ----echo = FALSE-------------------------------------------------------------
tbl <- tibble::tribble(
  ~Mathematics,   ~Physics,
  "function, $f (x)$",          "potential energy, $P(x)$",
  "derivative, $f\\prime (x)$",  "force, $F(x) = -P\\prime(x)$",
  "$\\min_x f(x) \\implies f\\prime (x) = 0$",
          "equilibrium, $\\min_x P(x) \\implies F (x) = 0$"
  )

tt(tbl)

## -----------------------------------------------------------------------------
# quadratic for potential energy

k <- 1/2
x <- seq(-6, 6, 0.5)
Px <- k * x^2

linepts <- function(a, b, x) {
  y <- a + b * x
  data.frame(x,y)
}

op <- par(mar = c(4,4,1,1)+.5)
plot(x, Px,
     type = "l",
     lwd = 2,
     axes = FALSE,
     ylab = "",
     cex.lab = 2)
axis(1)
axis(2)
abline(v=0, lwd = 1.2)
text(0, 18, "P(x)", cex = 2, xpd = TRUE)

# line through (3, 4.5), extending from x=2 to x=4
#   slope = d k X^2 / dx = X = 3
#   -> y - 4.5 = 3 * (x - 3) -> y = 3 X - 9
pts <- linepts(-4.5, 3, c(1.5, 3, 4.5))

lines(pts, col = "blue", lwd = 3)
points(pts[2,], pch = 16, cex= 2.3, col = "red")

text(-3, 14, expression(P(x) == frac(1,2) * X^2), cex = 1.2)
text(-2.8, 12, expression(F(x) == 2 * X), cex = 1.2)

points(0, 0, pch = 16, cex= 1.5, col = "red")
segments(x0 = -1.2, x1 = 1.2,
         y0 = 0, y1 = 0, col = "blue", lwd=2)

text(4 + c(0, .5), 6 + c(-0.7, -3),
     c(expression(F(x) == -dP/dx),
       expression(F(x) == - slope[x])),
     cex = 1.4,
     srt = 60,
     col = "blue",
     xpd = TRUE)
par(op)

## ----mean1--------------------------------------------------------------------
set.seed(1234)
N <- 8
df <- tibble(
  x = runif(N, 1, 10),
  y = seq(1, N)
)

means <- colMeans(df) 
xbar <- means[1] |> print()

## ----mean2--------------------------------------------------------------------
df <- df |>
  mutate(tension = abs(x - xbar),
         diameter = 0.2)

## ----mean-springs-------------------------------------------------------------
ggplot(df, aes(x=x, y=y)) +
  geom_point(size = 5, color = "red") +
  geom_segment(x = xbar, xend = xbar,
               y = 1/2,  yend = N + 1/2,
               linewidth = 3) +
  geom_spring(aes(x = x, xend = xbar,
                  y = y, yend = y,
                  tension = tension,
                  diameter = diameter),
              color = "blue",
              linewidth = 1.2) +
  labs(x = "Value (x)",
       y = "Observation number") +
  ylim(0, N+2) +
  scale_y_continuous(breaks = 1:N) +
  annotate("text", x = xbar, y = N + 1,
           label = "Movable\nrod", size = 5,
           lineheight = 3/4) 

## ----cent1--------------------------------------------------------------------
set.seed(1234)
N <- 10
df <- tibble(
  x = runif(N, 1, 10),
  y = runif(N, 1, 10)
)

## ----cent2--------------------------------------------------------------------
means <- colMeans(df)
xbar <- means[1]; ybar <- means[2]

## ----cent3--------------------------------------------------------------------
df <- df |>
  mutate(tension = sqrt((x - xbar)^2 + (y - ybar)^2),
         diameter = 0.4)

## ----centroid-springs---------------------------------------------------------
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
  scale_x_continuous(breaks = 1:10) +
  scale_y_continuous(breaks = 1:10) 

## -----------------------------------------------------------------------------
knitr::include_graphics(here::here("man/figures/levi-least-sq-springs.png"))

## -----------------------------------------------------------------------------
XPX <- matrix(c(
  "n", "\\sum x",
  "\\sum x", "\\sum x^2"), 2, 2)
XPY <- matrix(c("\\sum y", "\\sum x y"), 2,1)
b <- matrix(c("a", "b"), 2, 1)
Eqn(latexMatrix(XPX, matrix = "bmatrix"),
    latexMatrix(b),
    "& =",
    latexMatrix(XPY),
    Eqn_newline(),
    "\\mathbf{X}^\\top \\mathbf{X} \\mathbf{b}",
    ' & =',
    "\\mathbf{X}^\\top \\mathbf{y}",
    align = TRUE
)

## -----------------------------------------------------------------------------
knitr::include_graphics(here::here("man/figures/least-squares-geometry.png"))

## ----anscombe1----------------------------------------------------------------
simple_plot <-
  ggplot(anscombe,
         aes(x = x1, y = y1)) + 
  geom_point(size = 4, color = "red") +
  geom_smooth(method = lm, formula = y~x, 
              se = FALSE,
              linewidth = 1.4) 
simple_plot


## ----anscombe2----------------------------------------------------------------
model <- lm(y1 ~ x1, data = anscombe)

anscombe <- anscombe |>
  mutate(fitted = predict(model))

simple_plot +
  geom_spring(data=anscombe,
              aes(xend = x1, yend = fitted),
              color = "blue",
              linewidth = 1)

## ----augment------------------------------------------------------------------
broom::augment(model)

## ----anscombe3----------------------------------------------------------------
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

## -----------------------------------------------------------------------------
knitr::include_graphics(here::here("man/figures/Gaines-fig1.png"))

