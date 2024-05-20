#if the time t of psi_pa[t+1,] == Time, we evaluate psi_pa to be NA to make psi_tilde = 0
#psi_pa <- psi_pa[t+1,]
evaluate_psi_tilde <- function(x, psi_pa, model){  #from 0 to T. 0,T = 1 
  A <- model$A
  B <- model$B

  dif <- A%*%x - psi_pa[1:d]
    
  psi_tilde <- (-d/2)*log(2*pi) - (1/2)*log(det(diag(psi_pa[(d+1):(d+d)] + diag(B), nrow=d, ncol=d))) +
    (-1/2)*t(dif)%*%diag((psi_pa[(d+1):(d+d)] + diag(B))^(-1), nrow=d,ncol=d)%*%dif
  
  if(is.na(psi_tilde)){
    psi_tilde <- 0
  }
  
  return(psi_tilde)
}

#test:
#x <- c(0.7418182, 0.1584118)
#psi_pa <- c(-0.8921799, -0.4000272, 1.000000, 1.0000000)
#expected_output <- dmvn(as.vector(A%*%x), psi_pa[1:d], diag(psi_pa[(d+1):(d+d)] + diag(B), nrow=d, ncol=d), log = TRUE)
#result <- evaluate_psi_tilde(x, psi_pa, model)
#tolerance <- 1e-6
#is_close_enough <- abs(result - expected_output) < tolerance
#if (!is_close_enough) {
# stop("Test failed: Incorrect output from optimization function.")
#}
