#N=N[l], x = X_apf[t,,]

learn_psi <- function(x, obs, model){
  Time <- nrow(obs)
  N <- dim(x)[2]
  psi <- matrix(NA, nrow = Time, ncol = N)
  psi_pa <- matrix(NA, nrow = Time, ncol = 2*d)
  
  #calculate psi
  for(t in Time:1){
    print(t)
    if(t == Time){
      for(i in 1:N){
        dif <- x[t, i,] - obs[t,]
        
        psi[t,i] <- (1 / ((2 * pi)^(d / 2))) * 
          exp(-0.5 * t(dif) %*% dif)
      }
      
      
    }else{
      for(i in 1:N){
        
        psi[t,i] <- exp(evaluate_log_g(model, x[t,i,], obs[t,]))*
          exp(evaluate_psi_tilde(x[t,i,], psi_pa[t+1, ], model))
          
          #(2*pi)^(-d/2)*prod(psi_pa[t+1, (d+1):(d+d)]+1)^(-1/2)*
          #exp(-(1/2)*t(A%*%x[t,i,] - psi_pa[t+1, 1:d])%*%diag((psi_pa[t+1, (d+1):(d+d)]+diag(B))^(-1), nrow=d,ncol=d)%*%
           #     (A%*%x[t,i,] - psi_pa[t+1, 1:d]) )
        
      }
    }
    
    psi_pa[t,] <- optimization(x[t,,], log(psi[t,]))
    
    
    print(psi_pa[t, 1:d])
    print(obs[t,])
    
  }
  return(params = psi_pa)
  
}
