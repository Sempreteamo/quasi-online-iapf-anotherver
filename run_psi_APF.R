run_psi_APF <- function(n, psi_pa, N, L){ #purely filtering particles
  #l >= 2
  X <- array(NA, dim = c(Time, N, d))
  w <- matrix(NA, Time, N)
  logZ <- 0
  
  if(n == L){
    X[n-L+1,,] <- rnorm(N*d)  
    for (i in 1:(n-L)) {
      X[i,,] <- X[n-L+1,,]
    }
    for(i in 1:N){
      w_init[1:(n-L+1),i] <- g(obs[n-L+1,], X[n-L+1,i,])  
    }
  }else{
    #cat('n=',n,'L=',L)
    #cat('w[n-L]=', w[n-L,])
    output <- change_mu(w[n-L,], X[n-L,,])
    X[n-L+1,,] <- output[[1]]
    for (i in 1:(n-L)) {
      X[i,,] <- X[n-L+1,,]
    }
    w_ <- output[[2]]
    print('pass')
    
    
    for (i in 1:Num){
      # X[1:(n-L+1),i,] <- f(X[n-L,i,])
      w_init[1:(n-L+1), i] <- g(obs[n-L+1,], X[n-L+1,i,]) 
    }
    
  }
  
  for(t in (n-L+2):n){
    
    if(compute_ESS(w_init[t-1,], is.log=TRUE) <= kappa*N){
      
      ancestors <- resample(w[t-1,])
      logZ = logZ + normalise_weights_in_log_space(w[t-1,])[[2]]
      # at the initialization stage, we want filtering particles for psi
      
      for(i in 1:N){
        X[t,i,] <- f(X[t-1,mix[i],]) 
        w_init[t,i] <- g(obs[t,], X[t,i,])  
      }
      
    }else{
      for(i in 1:N){
        X[t,i,] <- f(X[t-1,i,]) 
        w_init[t,i] <- w_init[t-1,i] + g(obs[t,], X[t,i,])  
      }
    }
  }
  
  logZ <- logZ + normalise_weights_in_log_space(w[n,])[[2]]

  if(n == L){
    X[n-L+1,,] <- mu_aux(psi_pa, l, N, n-L+1)
    for (i in 1:(n-L)) {
      X[i,,] <- X[n-L+1,,]
    }
    for(i in 1:N){
      w[1:(n-L+1),i] <- g_aux(obs[n-L+1,], X[n-L+1,i,], n-L+1, psi_pa, n, L) 
    }
  }else{
    output <- change_mupsi(X[n-L,,], w[n-L,], psi_pa, n-L+1, N, l)
    X[n-L+1,,] <- output[[1]]
    for (i in 1:(n-L)) {
      X[i,,] <- X[n-L+1,,]
    }
    sum_ <- output[[2]]
    
    for (i in 1:N){
      #X[1:(n-L+1),i,] <- f_aux(X[n-L, i,], psi_pa, n-L+1)
      #w[1:(n-L+1), i] <- g_aux(obs[n-L+1,], X[n-L+1,i,], n-L+1, psi_pa, n, L)
      
      w[1:(n-L+1), i] <- g_transition(obs[n-L+1,], X[n-L+1,i,],  n-L+1, psi_pa, n) + log(sum_)
    }
    
  }
  
  for(t in (n-L+2):n){
    #print(t)
    if(compute_ESS(w[t-1,], is.log = TRUE) <= kappa*N){
      
      ancestors <- resample(w[t-1,])
      logZ = logZ + normalise_weights_in_log_space(w[t-1,])[[2]]
      
      for(i in 1:N){
        #filtering particles
        X[t,i,] <- f_aux(X[t-1, ancestors[i],], psi_pa, t)
        w[t,i] <- g_aux(obs[t,], X[t,i,], t, psi_pa, n, L) 
      }
    }else{
      
      for(i in 1:N){
        #filtering particles
        X[t,i,] <- f_aux(X[t-1,i,], psi_pa, t) 
        w[t,i] <- w[t-1,i] + g_aux(obs[t,], X[t,i,], t, psi_pa, n, L)
      }
    }
    
  }

  logZ <- logZ + normalise_weights_in_log_space(w[n,])[[2]]
  
  return(list(X, w, logZ, ancestors))
}
