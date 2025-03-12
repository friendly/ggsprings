# https://ggplot2-book.org/ext-springs

#' Create a spring
#'
#' @importFrom rlang abort
#' @param x      starting X coordinate
#' @param y      starting Y coordinate
#' @param xend   ending X coordinate
#' @param yend   ending Y coordinate
#' @param diameter diameter of spring
#' @param tension  tension of spring
#' @param n      number of points betwee start and end
#'
#' @return A data frame of x, y coordinates
#' @export
#'
#' @examples
#' spring <- create_spring(
#' x = 4, y = 2, xend = 10, yend = 6,
#' diameter = 2, tension = 0.6, n = 50
#' )
#'
#' ggplot(spring) +
#'   geom_path(aes(x = x, y = y)) +
#'   coord_equal()


create_spring <- function(x,
                          y,
                          xend,
                          yend,
                          diameter = 1,
                          tension = 0.75,
                          n = 50) {

  # Validate the input arguments
  if (tension <= 0) {
    rlang::abort("`tension` must be larger than zero.")
  }
  if (diameter == 0) {
    rlang::abort("`diameter` can not be zero.")
  }
  if (n == 0) {
    rlang::abort("`n` must be greater than zero.")
  }

  # Calculate the direct length of the spring path
  length <- sqrt((x - xend)^2 + (y - yend)^2)

  # Calculate the number of revolutions and points we need
  n_revolutions <- length / (diameter * tension)
  n_points <- n * n_revolutions

  # Calculate the sequence of radians and the x and y offset values
  radians <- seq(0, n_revolutions * 2 * pi, length.out = n_points)
  x <- seq(x, xend, length.out = n_points)
  y <- seq(y, yend, length.out = n_points)

  # Create and return the transformed data frame
  data.frame(
    x = cos(radians) * diameter/2 + x,
    y = sin(radians) * diameter/2 + y
  )
}

if(FALSE) {
  spring <- create_spring(
    x = 4, y = 2, xend = 10, yend = 6,
    diameter = 2, tension = 0.6, n = 50
  )

  ggplot(spring) +
    geom_path(aes(x = x, y = y)) +
    coord_equal()
}
