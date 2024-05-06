#params <- list(C, D)
evaluate_log_g <- function(params, x, datum){  
  C <- params[[1]]
  D <- params[[2]]
  dif <- datum - C%*%x
  log_density <- (-d/2)*log(2*pi) - (1/2)*log(det(D)) + (-1/2)*t(dif)%*%D%*%dif
  return (log_density)
}

#test
#x <- rnorm(d)
#datum <- rnorm(d)
