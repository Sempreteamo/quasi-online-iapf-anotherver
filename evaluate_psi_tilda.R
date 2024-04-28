evaluate_psi_tilda <- function(x, psi_pa, t, n){  #from 0 to T. 0,T = 1 
  if(t == n){
    psi_tilda <- 0
  }else{   #psi_pa_t = psi_t
    psi_tilda <- (-d/2)*log(2*pi) -(1/2)*log(det(diag(psi_pa[t+1, (d+1):(d+d)]+1, nrow=d, ncol=d))) +
      (-1/2)*t(A%*%x - psi_pa[t+1, 1:d])%*%diag((psi_pa[t+1, (d+1):(d+d)]+1)^(-1), nrow=d,ncol=d)%*%
      (A%*%x-psi_pa[t+1, 1:d]) 
  }
  return(psi_tilda)
}