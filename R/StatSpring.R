#' @name ggSpring
#'
#' @title ggSprings extensions to ggplot2
#' @format NULL
#' @usage NULL
#' @rdname Spring_protos
#' @export

StatSpring <- ggproto("StatSpring", Stat,

    setup_data = function(data, params) {
      if (anyDuplicated(data$group)) {
        data$group <- paste(data$group, seq_len(nrow(data)), sep = "-")
      }
      data
    },

    compute_panel = function(data, scales, n = 50) {
      cols_to_keep <- setdiff(names(data), c("x", "y", "xend", "yend"))
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
      do.call(rbind, springs)
    },

    required_aes = c("x", "y", "xend", "yend"),
    optional_aes = c("diameter", "tension")
)
