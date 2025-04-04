# Define StatSmoothFit
# from: https://github.com/EvaMaeRey/ggsmoothfit/blob/main/README.md#via-friendly-ggplot2-extenders-ggsprings-discussion-and-springs-extension-case-study

# Compute smooth predictions at the values of x observed in the data set or as specified by xseq=
# Also preserve the values of y (as yend) so that we can draw in the residual error.

compute_group_smooth_fit <- function(data, scales, method = NULL, formula = NULL,
                                     xseq = NULL,
                                     level = 0.95, method.args = list(),
                                     na.rm = FALSE, flipped_aes = NA){

  if(is.null(xseq)){ # predictions based on observations

    ggplot2::StatSmooth$compute_group(data = data, scales = scales,
                                      method = method, formula = formula,
                                      se = FALSE, n= 80, span = 0.75, fullrange = FALSE,
                                      xseq = data$x,
                                      level = .95, method.args = method.args,
                                      na.rm = na.rm, flipped_aes = flipped_aes) |>
      dplyr::mutate(xend = data$x,
                    yend = data$y)

  }else{  # predict specific input values

    ggplot2::StatSmooth$compute_group(data = data, scales = scales,
                                      method = method, formula = formula,
                                      se = FALSE, n= 80, span = 0.75, fullrange = FALSE,
                                      xseq = xseq,
                                      level = .95, method.args = method.args,
                                      na.rm = na.rm, flipped_aes = flipped_aes)

  }

}

# Pass to ggproto

StatSmoothFit <- ggplot2::ggproto("StatSmoothFit",
                                  ggplot2::StatSmooth,
                                  compute_group = compute_group_smooth_fit,
                                  required_aes = c("x", "y"))

# Pass to stat_*/ geom_ functions

#library(statexpress)

stat_smooth_fit <- function(geom = "point", ...){

  qlayer(geom = geom,
         stat = StatSmoothFit, ...)

}


geom_smooth_fit <- function(...){

  qlayer(geom = qproto_update(GeomPoint
                              # aes(colour = from_theme(accent))
                              ),
         stat = StatSmoothFit, ...)

}

geom_smooth_residuals <- function(...){

  qlayer(geom = qproto_update(GeomSegment
                              # aes(colour = from_theme(accent))
                              ),
         stat = StatSmoothFit, ...)

}

if (FALSE) {
mtcars %>%
  ggplot() +
  aes(x = wt, y = mpg) +
  geom_point() +
  geom_smooth() +
  geom_smooth_fit() +
  geom_smooth_residuals() +
  geom_smooth_fit(xseq = 2:3, size = 8)
}
