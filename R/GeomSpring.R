# https://ggplot2-book.org/ext-springs#creating-the-geom

#' @name ggSpring
#'
#' @title ggSprings extensions to ggplot2
#' @format NULL
#' @usage NULL
#' @importFrom rlang `%||%`
#' @rdname Spring_protos
#' @export

GeomSpring <- ggproto("GeomSpring", Geom,

    # Ensure that each row has a unique group id
    setup_data = function(data, params) {
      if (is.null(data$group)) {
        data$group <- seq_len(nrow(data))
      }
      if (anyDuplicated(data$group)) {
        data$group <- paste(data$group, seq_len(nrow(data)), sep = "-")
      }
      data
    },

    # Transform the data inside the draw_panel() method
    draw_panel = function(data,
                          panel_params,
                          coord,
                          n = 50,
                          arrow = NULL,
                          lineend = "butt",
                          linejoin = "round",
                          linemitre = 10,
                          na.rm = FALSE) {

      # Transform the input data to specify the spring paths
      cols_to_keep <- setdiff(names(data), c("x", "y", "xend", "yend"))

      # Set default for tension, diameter if not supplied
      data$diameter <- data$diameter %||% (.025 * abs(min(data$x) - max(data$x)))
      data$springlength <- sqrt((data$x - data$xend)^2 + (data$y - data$yend)^2)
      data$tension <-  data$tension %||% (1 * data$springlength)

      springs <- lapply(seq_len(nrow(data)), function(i) {
        spring_path <- create_spring(
          data$x[i],
          data$y[i],
          data$xend[i],
          data$yend[i],
          data$diameter[i],
          data$tension[i],
          n
        )
        cbind(spring_path, unclass(data[i, cols_to_keep]))
      })
      springs <- do.call(rbind, springs)

      # Use the draw_panel() method from GeomPath to do the drawing
      GeomPath$draw_panel(
        data = springs,
        panel_params = panel_params,
        coord = coord,
        arrow = arrow,
        lineend = lineend,
        linejoin = linejoin,
        linemitre = linemitre,
        na.rm = na.rm
      )
    },

    # Specify the default and required aesthetics
    required_aes = c("x", "y", "xend", "yend"),
    default_aes = aes(
      colour = "black",
      linewidth = 0.5,
      linetype = 1L,
      alpha = NA
    )
)
