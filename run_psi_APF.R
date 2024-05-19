run_psi_APF <- function(model, data, N, psi_pa, init = TRUE){ #purely filtering particles
  C <- model$C
  D <- model$D
  obs <- data[[1]]
  breaks <- data[[2]]
  Time <- nrow(obs)
    
    X <- array(NA, dim = c(Time, N, d))
    w <- matrix(NA, Time, N)
    logZ <- 0
    
    if(init){
      if(breaks[1] == 1){
        X[1,,] <- rnorm(N*d)  
        for(i in 1:N){
          w[1,i] <- evaluate_log_g(list(C, D), X[1,i,], obs[1,])  
        }
      }else{
        #cat('n=',n,'L=',L)
        #cat('w[n-L]=', w[n-L,])
        output <- change_mu(w[n-L,], X[n-L,,])
        X[1,,] <- output[[1]]
        w_ <- output[[2]]
        
        
        for (i in 1:Num){
          w[1, i] <- evaluate_log_g(list(C, D), X[1,i,], obs[1,])
        }
        
      }
      
      for(t in 2:Time){
        
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
      if(breaks[1] == 1){
        X[1,,] <- sample_twisted_initial(list(model$ini_mu, model$ini_cov), psi_pa, N)
        for(i in 1:N){
          w[1,i] <- eval_twisted_potential(model, psi_pa, X[1,i,], obs[1,])
        }
      }else{
        output <- change_mupsi(X[n-L,,], w[n-L,], psi_pa, 1, N, l)
        X[1,,] <- output[[1]]
        
        for (i in 1:N){
          #X[1:(1),i,] <- f_aux(X[n-L, i,], psi_pa, 1)
          #w[1:(1), i] <- g_aux(obs[1,], X[1,i,], 1, psi_pa, n, L)
          
          w[1, i] <- eval_twisted_potential(model, psi_pa, X[1,i,], obs[1,])
        }
        
      }
      
      for(t in 2:Time){
        #print(t)
        if(compute_ESS_log(w[t-1,]) <= kappa*N){
          
          ancestors <- resample(w[t-1,])
          logZ = logZ + normalise_weights_in_log_space(w[t-1,])[[2]]
          
          for(i in 1:N){
            #filtering particles
            X[t,i,] <- sample_twisted_transition(X[t-1, ancestors[i],], list(model$A, model$B), psi_pa, N)
            w[t,i] <- eval_twisted_potential(model, psi_pa, X[t,i,], obs[t,])
            
          }
        }else{
          
          for(i in 1:N){
            #filtering particles
            X[t,i,] <- sample_twisted_transition(X[t-1, i,], list(model$A, model$B), psi_pa, N)
            w[t,i] <- w[t-1,i] + eval_twisted_potential(model, psi_pa, X[t,i,], obs[t,])
          }
        }
        
      }
      
      logZ <- logZ + normalise_weights_in_log_space(w[n,])[[2]] 
    }
    
    
  
  
  return(list(X, w, logZ, ancestors))
}
