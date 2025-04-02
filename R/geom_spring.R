# constructors
# https://ggplot2-book.org/ext-springs#a-constructor

#' Connect observations with springs
#'
#' \code{geom_spring} is similar to \code{\link[ggplot2]{geom_path}} in that it connects points,
#' but uses a spring instead of a line.
#'
#' @inheritParams ggplot2::geom_path
#' @param n Number of points
#'
#' @section Aesthetics:
#' geom_spring understands the following aesthetics (required aesthetics are in bold):
#'
#' - **x**
#' - **y**
#' - **xend**
#' - **yend**
#' - diameter
#' - tension
#' - color
#' - linewidth
#' - linetype
#' - alpha
#' - lineend
#'
#' The additional aesthetics are explained below:
#' \describe{
#'   \item{diameter}{Diameter of the spring, i.e., the diameter of a circle
#'   that is stretched into a spring shape.}
#'   \item{tension}{Spring tension constant. This is calibrated as the total
#'   distance moved from the start point to the end point, divided by the size
#'   of the generating circle.}
#' }
#'
#' @importFrom ggplot2 layer
#'
#' @return A ggplot2 layer
#' @export
#'
#' @examples
#' # None yet
geom_spring <- function(mapping = NULL,
                        data = NULL,
                        stat = "identity",
                        position = "identity",
                        ...,
                        n = 50,
                        arrow = NULL,
                        lineend = "butt",
                        linejoin = "round",
                        na.rm = FALSE,
                        show.legend = NA,
                        inherit.aes = TRUE) {
  layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomSpring,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      arrow = arrow,
      lineend = lineend,
      linejoin = linejoin,
      na.rm = na.rm,
      ...
    )
  )
}

#' Spring stat
#'
#' @inheritParams ggplot2::geom_path
#' @importFrom ggplot2 layer
#  @param mapping
#  @param data
#' @param geom     The \code{geom} used to draw the spring segment
#  @param position
#  @param ...
#' @param diameter Diameter of the spring, i.e., the diameter of a circle that is stretched into a spring shape.
#' @param tension  Spring tension constant. This is calibrated as the total distance moved from the start point to the end point, divided by the size of the generating circle.
#' @param n        Number of points
#  @param na.rm
#  @param show.legend
#  @param inherit.aes
#'
#' @return A ggplot2 layer
#' @export
#'
#' @examples
#' # None yet
stat_spring <- function(mapping = NULL,
                        data = NULL,
                        geom = "path",
                        position = "identity",
                        ...,
                        diameter = 1,
                        tension = 0.75,
                        n = 50,
                        na.rm = FALSE,
                        show.legend = NA,
                        inherit.aes = TRUE) {
  layer(
    data = data,
    mapping = mapping,
    stat = StatSpring,
    geom = geom,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      diameter = diameter,
      tension = tension,
      n = n,
      na.rm = na.rm,
      ...
    )
  )
}
