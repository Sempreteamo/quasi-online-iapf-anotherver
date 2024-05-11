#if the time n == the time t of psi_pa[t,] we evaluate psi_pa to be zero? to make psi_tilda = 0
#psi_pa <- psi_pa[t,]
evaluate_psi_tilda <- function(x, psi_pa, model){  #from 0 to T. 0,T = 1 
  A <- model$A
  B <- model$B

  dif <- A%*%x - psi_pa[1:d]
    
  psi_tilda <- (-d/2)*log(2*pi) - (1/2)*log(det(diag(psi_pa[(d+1):(d+d)] + diag(B), nrow=d, ncol=d))) +
    (-1/2)*t(dif)%*%diag((psi_pa[(d+1):(d+d)] + diag(B))^(-1), nrow=d,ncol=d)%*%dif
  
  return(psi_tilda)
}

#test:psi_pa <- rnorm(d*2)
#x <- rnorm(d)
#evaluate_psi_tilda(x, psi_pa, model)
