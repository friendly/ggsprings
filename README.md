
# ggsprings

<!-- badges: start -->
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

ggsprings is designed to
implement an extension of `geom_path` which draws paths as springs instead of straight lines.
Aside from possible artistic use, the main impetus for this is to draw points connected by springs,
with properties of length and tension. 

A leading example is to illustrate how least squares regression
is "solved" by connecting data points to a rod, where the springs are constrained to be vertical.
The mathematics behind this are well-described in this [Math Stackexchange post](https://math.stackexchange.com/questions/2369673/proving-linear-regression-by-using-physical-springs-model),
where the least squares estimates of intercept and slope are shown to be the equilibrium position that minimized the sum of forces
and torques exerted by springs.

![](man/figures/potential-energy.png)

If the springs are allowed to be free, the physical solution is the major PCA axis.


How to do this is described in the `ggplot2` book, https://ggplot2-book.org/ext-springs.
The current (non-working) version here was copied/pasted from the book.

A blog post by Joshua Loftus, [Least squares by springs](https://joshualoftus.com/posts/2020-11-23-least-squares-as-springs/least-squares-as-springs.html)
illustrates this, citing [code from Thomas Lin Pederson](https://twitter.com/thomasp85/status/1331338379636649986).
Code to reproduce the first example is contained in `examples/springs.R` and `examples/gapminder-ex.R`.

### Examples

These images show the intent of a `ggsprings` package.

**Least squares regression**

A plot of `lifeExp` vs. `gdpPercap` from the `gapminder` data, with `gdpPercap` on a log10 scale, using the code in the `examples/` folder.
Springs are connected between the observed value `y = lifeExp` and the fitted value on the regression line, `yend = yhat`, computed
with `predict()` for the linear model.
`tension` was set to `5 + (lifeExp - yhat)^2)`

```
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
```

![](man/figures/loftus-springs-ex1.png)

**Principal components analysis**

![](man/figures/loftus-springs-ex2.png)

**Animated version**

This [StatsExchange post](https://stats.stackexchange.com/questions/2691/making-sense-of-principal-component-analysis-eigenvectors-eigenvalues/140579#140579)
show an animation of the process of fitting PCA by springs.

![](man/figures/pca-springs-cropped.gif)


## Installation

You can install the development version of ggsprings by forking this repo and trying to get it going.


## Example

There is really nothing working here yet:

``` r
library(ggsprings)
## basic example code
```

## Related 

* An [interactive demo](https://www.desmos.com/calculator/90vaqtqpx6) by Trey Goesh allows you to 
visualize the effect of moving points, changing spring parameters, etc.
