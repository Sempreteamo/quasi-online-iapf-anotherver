evaluate_psi_tilda <- function(x, psi_pa, t, n){  #from 0 to T. 0,T = 1 
  if(t == n){
    psi_tilda <- 0
  }else{   #psi_pa_t = psi_t
    dif <- A%*%x - psi_pa[t+1, 1:d]
    
    psi_tilda <- (-d/2)*log(2*pi) -(1/2)*log(det(diag(psi_pa[t+1, (d+1):(d+d)]+1, nrow=d, ncol=d))) +
      (-1/2)*t(dif)%*%diag((psi_pa[t+1, (d+1):(d+d)]+1)^(-1), nrow=d,ncol=d)%*%dif
  }
  return(psi_tilda)
}
