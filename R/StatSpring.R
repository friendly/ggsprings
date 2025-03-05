StatSpring <- ggproto("StatSpring", Stat,

    # Edit the input data to ensure the group identifiers are unique
    setup_data = function(data, params) {
      if (anyDuplicated(data$group)) {
        data$group <- paste(data$group, seq_len(nrow(data)), sep = "-")
      }
      data
    },

    # Construct data for this panel by calling create_spring()
    compute_panel = function(data,
                             scales,
                             diameter = 1,
                             tension = 0.75,
                             n = 50) {
      cols_to_keep <- setdiff(names(data), c("x", "y", "xend", "yend"))
      springs <- lapply(
        seq_len(nrow(data)),
        function(i) {
          spring_path <- create_spring(
            data$x[i],
            data$y[i],
            data$xend[i],
            data$yend[i],
            diameter = diameter,
            tension = tension,
            n = n
          )
          cbind(spring_path, unclass(data[i, cols_to_keep]))
        }
      )
      do.call(rbind, springs)
    },

    # Specify which aesthetics are required input
    required_aes = c("x", "y", "xend", "yend")
)
