#params <- list(C, D)
evaluate_log_g <- function(model, x, datum){  
  C <- model$C
  D <- model$D
  dif <- datum - C%*%x
  log_density <- (-d/2)*log(2*pi) - (1/2)*log(prod(diag(D))) + (-1/2)*t(dif)%*%D%*%dif
  return (log_density) 
}

#test
#x <- c(0.7418182, 0.1584118)
#datum <- c(0.67968786, -0.01167335)
#expected_output <- dmvn(datum, C%*%x, D, log = TRUE)
#result <- evaluate_log_g(model, x, datum)
#tolerance <- 1e-6
#is_close_enough <- abs(result - expected_output) < tolerance
#if (!is_close_enough) {
# stop("Test failed: Incorrect output from optimization function.")
#}
