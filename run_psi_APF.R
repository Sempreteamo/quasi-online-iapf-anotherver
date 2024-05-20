run_psi_APF <- function(model, data, N, psi_pa, init){ #purely filtering particles
  A <- model$A
  B <- model$B
  C <- model$C
  D <- model$D
  obs <- data[[1]]
  breaks <- data[[2]]
  w_previous <- data[[3]]
  X_previous <- data[[4]]
  Time <- nrow(obs)
    
    X <- array(NA, dim = c(Time, N, d))
    w <- matrix(NA, Time, N)
    logZ <- 0
    
    if(init){
      if(breaks[1] == 1){
        X[1,,] <- rnorm(N*d)  
        for(i in 1:N){
          w[1,i] <- evaluate_log_g(model, X[1,i,], obs[1,])  
        }
      }else{
        
        s <- resample(w_previous, mode = 'multi')
        for(i in 1:N){
          X[1,i,] <- rmvn(1, A%*%X_previous[s[i],], B)
        }
        
        for (i in 1:N){
          w[1, i] <- evaluate_log_g(model, X[1,i,], obs[1,])
        }
        
      }
      
      for(t in 2:Time){
        
        if(compute_ESS_log(w[t-1,]) <= kappa*N){
          
          ancestors <- resample(w[t-1,])
          logZ = logZ + normalise_weights_in_log_space(w[t-1,])[[2]]
          # at the initialization stage, we want filtering particles for psi
          
          X[t,,] <- sample_normal_distribution(list(t(A%*%t(X[t-1, ancestors,])), B), N)
          for(i in 1:N){
            w[t,i] <- evaluate_log_g(model, X[t,i,], obs[t,])  
          }
          
        }else{
          
          X[t,,] <- sample_normal_distribution(list(t(A%*%t(X[t-1, ,])), B), N)
          for(i in 1:N){
            w[t,i] <- w[t-1,i] + evaluate_log_g(model, X[t,i,], obs[t,])  
          }
        }
      }
      
      logZ <- logZ + normalise_weights_in_log_space(w[Time,])[[2]]
      
    }else{
      if(breaks[1] == 1){
        X[1,,] <- sample_twisted_initial(list(model$ini_mu, model$ini_cov), psi_pa[1,], N)
        for(i in 1:N){
          w[1,i] <- eval_twisted_potential(model, list(NA, psi_pa[2,], psi_pa[1,]), X[1,i,], obs[1,])
        }
      }else{
        output <- change_mupsi(X[n-L,,], w_previous, psi_pa, 1, N, l)
        X[1,,] <- output[[1]]
        
        for (i in 1:N){
          #X[1:(1),i,] <- f_aux(X[n-L, i,], psi_pa, 1)
          #w[1:(1), i] <- g_aux(obs[1,], X[1,i,], 1, psi_pa, n, L)
          
          w[1, i] <- eval_twisted_potential(model, list(NA, psi_pa[2,], psi_pa[1,]), X[1,i,], obs[1,])
        }
        
      }
      
      for(t in 2:(Time-1)){
        #print(t)
        if(compute_ESS_log(w[t-1,]) <= kappa*N){
          
          ancestors <- resample(w[t-1,])
          logZ = logZ + normalise_weights_in_log_space(w[t-1,])[[2]]
          
          for(i in 1:N){
            #filtering particles
            X[t,i,] <- sample_twisted_transition(X[t-1, ancestors[i],], list(model$A, model$B), psi_pa[t,], 1)
            w[t,i] <- eval_twisted_potential(model, list(psi_pa[t,], psi_pa[t+1,], psi_pa[t,]), X[t,i,], obs[t,])
            
          }
        }else{
          
          for(i in 1:N){
            #filtering particles
            X[t,i,] <- sample_twisted_transition(X[t-1, i,], list(model$A, model$B), psi_pa[t,], 1)
            w[t,i] <- w[t-1,i] + eval_twisted_potential(model, list(psi_pa[t,], psi_pa[t+1,], psi_pa[t,]), X[t,i,], obs[t,])
          }
        }
        
      }
      
      t = Time
      
      if(compute_ESS_log(w[t-1,]) <= kappa*N){
        
        ancestors <- resample(w[t-1,])
        logZ = logZ + normalise_weights_in_log_space(w[t-1,])[[2]]
        
        for(i in 1:N){
          #filtering particles
          X[t,i,] <- sample_twisted_transition(X[t-1, ancestors[i],], list(model$A, model$B), psi_pa[t,], 1)
          w[t,i] <- eval_twisted_potential(model, list(psi_pa[t,], NA, psi_pa[t,]), X[t,i,], obs[t,])
          
        }
      }else{
        
        for(i in 1:N){
          #filtering particles
          X[t,i,] <- sample_twisted_transition(X[t-1, i,], list(model$A, model$B), psi_pa[t,], 1)
          w[t,i] <- w[t-1,i] + eval_twisted_potential(model, list(psi_pa[t,], NA, psi_pa[t,]), X[t,i,], obs[t,])
        }
      }
      
      logZ <- logZ + normalise_weights_in_log_space(w[Time,])[[2]] 
    }
    
    
  
  
  return(list(X, w, logZ, ancestors))
}
