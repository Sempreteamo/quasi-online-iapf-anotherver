#if when the t of psi_pa[t,] = Time+1, where Time is the n we're evaluating at, then psi_pa == 0 to make 
#psi_t <- 0
#psi_pa <- psi_pa[t,]
evaluate_psi <- function(x, psi_pa){ 
  
    dif <- x - psi_pa[1:d]
    psi_t <- -(d/2)*log(2*pi) - (1/2)*log(det(diag(psi_pa[(d+1):(d+d)], nrow=d,ncol=d))) +						
      (-1/2)*t(dif)%*%diag((psi_pa[(d+1):(d+d)])^(-1), nrow=d,ncol=d)%*%dif
  
  return(psi_t)
}
