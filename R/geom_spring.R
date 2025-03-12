# constructors
# https://ggplot2-book.org/ext-springs#a-constructor

#' Connect observations with springs
#'
#' \code{geom_spring} is similar to \code{\link[ggplot2]{geom_path}} in that it connects points,
#' but uses a spring instead of a line.
#'
#' @inheritParams ggplot2::geom_path
#' @importFrom ggplot2 layer
#  @param mapping
#  @param data
#  @param stat
#  @param position
#  @param ...
#' @param diameter Diameter of the spring
#' @param tension  Spring tension constant
#' @param n        Number of points
#  @param arrow
#  @param lineend
#  @param linejoin
#  @param na.rm
#  @param show.legend
#  @param inherit.aes
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
      n = n,
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
#  @param geom
#  @param position
#  @param ...
#' @param n        Number of points
#  @param na.rm
#  @param show.legend
#  @param inherit.aes
#'
#' @return A ggplot2 layer
#' @export
#'
#' @examples
#' None yet
stat_spring <- function(mapping = NULL,
                        data = NULL,
                        geom = "path",
                        position = "identity",
                        ...,
                        # diameter = 1,
                        # tension = 0.75,
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
