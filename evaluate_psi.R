evaluate_psi <- function(x, psi_pa, Time){ 
  if(t == (Time + 1)){
    psi_t <- 0
  }else{
    dif <- x - psi_pa[t, 1:d]
    psi_t <- -(d/2)*log(2*pi) - (1/2)*log(det(diag(psi_pa[t, (d+1):(d+d)], nrow=d,ncol=d))) +						
      (-1/2)*t(dif)%*%diag((psi_pa[t, (d+1):(d+d)])^(-1), nrow=d,ncol=d)%*%dif
  }
  return(psi_t)
}
