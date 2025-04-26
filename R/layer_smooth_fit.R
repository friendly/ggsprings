StatSmoothFit <- ggplot2::ggproto("StatSmoothFit", 
                                  ggplot2::StatSmooth,
                                  compute_group = compute_group_smooth_fit,
                                  required_aes = c("x", "y"))

aes_color_accent <- GeomSmooth$default_aes[c("colour")]

GeomPointAccent <- ggproto("GeomPointAccent", GeomPoint, 
              default_aes = modifyList(GeomPoint$default_aes, 
                                       aes_color_accent))

GeomSegmentAccent <- ggproto("GeomSegmentAccent", GeomSegment,
                           default_aes = modifyList(GeomSegment$default_aes, 
                                                    aes_color_accent))

GeomSpringAccent <- ggproto("GeomSpringAccent", GeomSpring,
                           default_aes = modifyList(GeomSpring$default_aes,
                                                    aes_color_accent))


#' @export
layer_smooth_fit <- function (mapping = NULL, data = NULL, stat = StatSmoothFit, geom = GeomPointAccent, position = "identity", 
    ..., show.legend = NA, inherit.aes = TRUE) 
{
    layer(data = data, mapping = mapping, stat = stat, 
        geom = geom, position = position, show.legend = show.legend, 
        inherit.aes = inherit.aes, params = rlang::list2(na.rm = FALSE, 
            ...))
}

#' @export
stat_smooth_fit <- function(...){layer_smooth_fit(stat = StatSmoothFit, ...)}

#' @export
geom_smooth_fit <- function(...){layer_smooth_fit(geom = GeomPointAccent, ...)}

#' @export
geom_residuals <- function(...){layer_smooth_fit(geom = GeomSegmentAccent, ...)}

#' @export
geom_residual_springs <- function(...){layer_smooth_fit(geom = GeomSpringAccent, ...)}

