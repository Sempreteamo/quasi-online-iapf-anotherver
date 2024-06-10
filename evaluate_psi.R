#if when the t of psi_pa[t,] = Time+1, where Time is the n we're evaluating at, then psi_pa == NA to make 
#psi_x <- 0
#psi_pa <- psi_pa[t,]
evaluate_psi <- function(x, psi_pa){ 
    d <- length(x)
    dif <- x - psi_pa[1:d]
    psi_x <- -(d/2)*log(2*pi) - (1/2)*log(det(diag(psi_pa[(d+1):(d+d)], nrow=d,ncol=d))) +						
      (-1/2)*t(dif)%*%diag((psi_pa[(d+1):(d+d)])^(-1), nrow=d,ncol=d)%*%dif
  
  if(is.na(psi_x)){
    psi_x <- 0
  }
  
    return(psi_x) 
  
}
#test:
#x <- c(0.7418182, 0.1584118)
#psi_pa <- c(-0.8921799, -0.4000272, 1.000000, 1.0000000)
#expected_output <- dmvn(x, psi_pa[1:d], diag(psi_pa[(d+1):(d+d)], d, d), log = TRUE)
#result <- evaluate_psi(x, psi_pa)
#tolerance <- 1e-6
#is_close_enough <- abs(result - expected_output) < tolerance
#if (!is_close_enough) {
 # stop("Test failed: Incorrect output from optimization function.")
#}
