
# ggsprings

<!-- badges: start -->
<!-- badges: end -->

ggsprings is designed to
implement an extension of `geom_path` which draws paths as springs instead of straight lines.
Aside from possible artistic use, the main impetus for this is to draw points connected by springs,
with properties of length and tension. 

A leading example is to illustrate how least squares regression
is "solved" by connecting data points to a rod, where the springs are constrained to be vertical.
If the springs are allowed to be free, the physical solution is the major PCA axis.

How to do this is described in the `ggplot2` book, https://ggplot2-book.org/ext-springs.
The current (non-working) version here was copied/pasted from the book.

A blog post by Joshua Loftus, [Least squares by springs](https://ggplot2-book.org/ext-springs)
illustrates this, citing [code from Thomas Lin Pederson](https://twitter.com/thomasp85/status/1331338379636649986).

## Installation

You can install the development version of ggsprings by forking this repo and trying to get it going.


## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(ggsprings)
## basic example code
```

