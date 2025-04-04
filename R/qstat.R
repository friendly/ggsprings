#' @export
qstat <- function(compute_group = ggplot2::Stat$compute_group, ...) {

  ggplot2::ggproto("StatTemp", Stat, compute_group = compute_group, ...)

}

#' @export
qstat_group <- function(qstat_group, ...) {

  ggplot2::ggproto("StatTemp", Stat, qstat_group = qstat_group, ...)

}

#' @export
qstat_panel <- function(compute_panel, ...) {

  ggplot2::ggproto("StatTemp", Stat, compute_panel = compute_panel, ...)

}

#' @export
qstat_layer <- function(compute_layer, ...) {

  ggplot2::ggproto("StatTemp", Stat, compute_layer = compute_layer, ...)

}


#' title Quickly add a layer to a ggplot
#'
#' @export
qlayer <- function (mapping = NULL,
                    data = NULL,
                    geom = GeomPoint,
                    stat = StatIdentity,
                    position = position_identity(),
                    ...,
                    na.rm = FALSE,
                    show.legend = NA,
                    inherit.aes = TRUE)
{
  ggplot2::layer(
    data = data,
    mapping = mapping,
    geom = geom,
    stat = stat,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = rlang::list2(na.rm = na.rm, ...)
  )
}

#' @export
proto_update <- function(`_class`, `_inherit`, default_aes_update = NULL, ...){

  if(!is.null(default_aes_update)){

    default_aes <- aes(!!!modifyList(`_inherit`$default_aes, default_aes_update))

  }

  ggplot2::ggproto(`_class` = `_class`,
                   `_inherit` = `_inherit`,
                   default_aes = default_aes, ...)

}

#' @export
qproto_update <- function(`_inherit`, default_aes_update = NULL, ...){

  proto_update("protoTemp",
               `_inherit`,
               default_aes_update = default_aes_update,
               ...)
}
