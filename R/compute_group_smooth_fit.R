compute_group_smooth_fit <- function(data, scales, method = NULL, formula = NULL,
                           xseq = NULL,
                           level = 0.95, method.args = list(),
                           na.rm = FALSE, flipped_aes = NA){
  
  if(is.null(xseq)){ # predictions based on observations 

  StatSmooth$compute_group(data = data, scales = scales, 
                       method = method, formula = formula, 
                       se = FALSE, n= 80, span = 0.75, fullrange = FALSE,
                       xseq = data$x, 
                       level = .95, method.args = method.args, 
                       na.rm = na.rm, flipped_aes = flipped_aes) |>
      dplyr::mutate(xend = data$x,
                    yend = data$y)
  
  }else{  # predict specific input values
    
  StatSmooth$compute_group(data = data, scales = scales, 
                       method = method, formula = formula, 
                       se = FALSE, n= 80, span = 0.75, fullrange = FALSE,
                       xseq = xseq, 
                       level = .95, method.args = method.args, 
                       na.rm = na.rm, flipped_aes = flipped_aes)   
    
  }
  
}
