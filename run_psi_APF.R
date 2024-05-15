run_psi_APF <- function(model, data, N, psi_pa, init = TRUE, block, index){ #purely filtering particles
  obs <- data$obs
  breaks <- data$breaks
  psi_index <- data$psi_index
  C <- model$C
  D <- model$D
  Time <- nrow(obs)
    
    X <- array(NA, dim = c(Time, N, d))
    w <- matrix(NA, Time, N)
    logZ <- 0
    
    if(init){
      if(block){
        X[1,,] <- rnorm(N*d)  
        for(i in 1:N){
          w[1,i] <- evaluate_log_g(list(C, D), X[1,i,], obs[1,])  
        }
      }else{
        #cat('n=',n,'L=',L)
        #cat('w[n-L]=', w[n-L,])
        output <- change_mu(w[n-L,], X[n-L,,])
        X[1,,] <- output[[1]]
        for (i in 1:(n-L)) {
          X[i,,] <- X[1,,]
        }
        w_ <- output[[2]]
        
        
        for (i in 1:Num){
          # X[1:(1),i,] <- f(X[n-L,i,])
          w[1, i] <- g(obs[1,], X[1,i,]) 
        }
        
      }
      
      for(t in 2:(Time)){
        
        if(compute_ESS_log(w[t-1,]) <= kappa*N){
          
          ancestors <- resample(w[t-1,])
          logZ = logZ + normalise_weights_in_log_space(w[t-1,])[[2]]
          # at the initialization stage, we want filtering particles for psi
          
          X[t,,] <- sample_normal_distribution(list(t(A%*%t(X[t-1, ancestors,])), B), N)
          for(i in 1:N){
            w[t,i] <- evaluate_log_g(list(C, D), obs[t,], X[t,i,])  
          }
          
        }else{
          
          X[t,,] <- sample_normal_distribution(list(t(A%*%t(X[t-1, ,])), B), N)
          for(i in 1:N){
            w[t,i] <- w[t-1,i] + evaluate_log_g(list(C, D), obs[t,], X[t,i,])  
          }
        }
      }
      
      logZ <- logZ + normalise_weights_in_log_space(w[Time,])[[2]]
      
    }else{
      if(n == L){
        X[1,,] <- mu_aux(psi_pa, l, N, 1)
        for (i in 1:(n-L)) {
          X[i,,] <- X[1,,]
        }
        for(i in 1:N){
          w[1:(1),i] <- g_aux(obs[1,], X[1,i,], 1, psi_pa, n, L) 
        }
      }else{
        output <- change_mupsi(X[n-L,,], w[n-L,], psi_pa, 1, N, l)
        X[1,,] <- output[[1]]
        for (i in 1:(n-L)) {
          X[i,,] <- X[1,,]
        }
        sum_ <- output[[2]]
        
        for (i in 1:N){
          #X[1:(1),i,] <- f_aux(X[n-L, i,], psi_pa, 1)
          #w[1:(1), i] <- g_aux(obs[1,], X[1,i,], 1, psi_pa, n, L)
          
          w[1:(1), i] <- g_transition(obs[1,], X[1,i,],  1, psi_pa, n) + log(sum_)
        }
        
      }
      
      for(t in (n-L+2):n){
        #print(t)
        if(compute_ESS_log(w[t-1,], is.log = TRUE) <= kappa*N){
          
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
    }
    
    
  
  
  return(list(X, w, logZ, ancestors))
}

