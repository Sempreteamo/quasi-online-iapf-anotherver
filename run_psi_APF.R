#' Function of psi_APF
#' 
#' This function performs sampling and resampling procedures of 
#' psi-Auxiliary Particle Filter with κ-adaptive resampling
#'
#' @param model List containing model parameters
#' @param data List containing time series on which to run the filter
#' @param N Number of particles to use
#' @param psi_pa Parameters of the twisting function
#' @param init 
#'
#' @return A list containing:
#' x is the particle set generated during the sampling and resampling approach
#' w is the weights of particles
#' ancestors is indices associated with every resampling event
#' logZ is normalizing constant estimate
#' 
#' @export
#'
run_psi_APF <- function(model, data, N, psi_pa, init){ 
  ini_mu <- model$ini_mu
  ini_cov <- model$ini_cov
  A <- model$A
  B <- model$B
  C <- model$C
  D <- model$D
  d <- model$d
  obs <- data[[1]]
  breaks <- data[[2]]
  w_previous <- data[[3]]
  X_previous <- data[[4]]
  Time <- nrow(obs)
  kappa <- model$kappa
  ancestors <- matrix(NA, Time, N)
  
  X <- array(NA, dim = c(Time, N, d))
  w <- matrix(NA, Time, N)
  logZ <- 0
  
  #'init' controls the the condition where we conduct the initialization
  if(init){
    if(breaks[1] == 1){
      #the first block. break controls which block the algorithm is running
      
      X[1,,] <- rnorm(N*d)  
      for(i in 1:N){
        w[1,i] <- evaluate_log_g(model, X[1,i,], obs[1,])  
      }
      
      
    }else{
      
      s <- resample(w_previous, mode = 'multi')
      for(i in 1:N){
        X[1,i,] <- rmvn(1, A%*%X_previous[s[i],], B)
        w[1, i] <- evaluate_log_g(model, X[1,i,], obs[1,])
      }
    }
    
    ancestors[1,] <- seq(1:N) 
    
    for(t in 2:Time){
      if(compute_ESS_log(w[t-1,]) <= kappa*N){
        
        ancestors[t,] <- resample(w[t-1,], mode = 'multi')
        logZ = logZ + normalise_weights_in_log_space(w[t-1,])[[2]]
        
        for(i in 1:N){
          
          X[t,i,] <- rmvn(1, A%*%X[t-1, ancestors[t,i],], B)
          w[t,i] <- evaluate_log_g(model, X[t,i,], obs[t,])  
        }
        
      }else{
        ancestors[t,] <- ancestors[t-1,] 
        for(i in 1:N){
          X[t,i,] <- rmvn(1, A%*%X[t-1, i,], B)
          w[t,i] <- w[t-1,i] + evaluate_log_g(model, X[t,i,], obs[t,])  
        }
      }
    }
    
    
    logZ <- logZ + normalise_weights_in_log_space(w[Time,])[[2]]
    
  }else{
    if(breaks[1] == 1){
    
      X[1,,] <- sample_twisted_initial(list(mean = ini_mu, cov = ini_cov), psi_pa[1,], N)
      
      for(i in 1:N){
        w[1,i] <- eval_twisted_potential(model, list(psi_pa[1,], psi_pa[2,], psi_pa[1,]), X[1,i,], obs[1,])
      }
      
    }else{
      
      w_adj <- vector()
      
      for(i in 1:N){
        w_adj[i] <- w_previous[i]*exp(1/2*(t(A%*%X_previous[i,]) + t(psi_pa[1,1:d])%*%
                                             diag(psi_pa[1, (d+1):(d+d)]^(-2), nrow=d,ncol=d))%*%diag((psi_pa[1, (d+1):(d+d)]^(-2) + 1)^(-1), nrow=d,ncol=d)%*%
                                        (A%*%X_previous[i,] + diag(psi_pa[1, (d+1):(d+d)]^(-2), nrow=d,ncol=d)%*%psi_pa[1,1:d]) - 1/2*(t(A%*%X_previous[i,])%*%A%*%X_previous[i,] +
                                                                                                                                         t(psi_pa[1,1:d])%*%diag(psi_pa[1, (d+1):(d+d)]^(-2), nrow=d,ncol=d)%*%psi_pa[1,1:d]))
      }
      
      s <- resample(w_adj, mode = 'multi')
      
      for (i in 1:N){
        X[1, i, ] <- sample_twisted_transition(X_previous[s[i],], model, psi_pa[1,], 1)
        w[1, i] <- eval_twisted_potential(model, list(psi_pa[1,], psi_pa[2,], psi_pa[1,]), X[1,i,], obs[1,])
      }
      
    }
    
    ancestors[1,] <- seq(1:N) 
    
    for(t in 2:(Time-1)){
      if(compute_ESS_log(w[t-1,]) <= kappa*N){
        
        ancestors[t,] <- resample(w[t-1,], mode = 'multi')
        logZ = logZ + normalise_weights_in_log_space(w[t-1,])[[2]]
        
        for(i in 1:N){
          
          X[t,i,] <- sample_twisted_transition(X[t-1, ancestors[t,i],], model, psi_pa[t,], 1)
          w[t,i] <- eval_twisted_potential(model, list(NA, psi_pa[t+1,], psi_pa[t,]), X[t,i,], obs[t,])
          
        }
      }else{
        ancestors[t,] <- ancestors[t-1,]
        for(i in 1:N){
          
          X[t,i,] <- sample_twisted_transition(X[t-1, i,], model, psi_pa[t,], 1)
          w[t,i] <- w[t-1,i] + eval_twisted_potential(model, list(NA, psi_pa[t+1,], psi_pa[t,]), X[t,i,], obs[t,])
        }
      }
      
    }
   
    t = Time
    
    if(compute_ESS_log(w[t-1,]) <= kappa*N){
      
      ancestors[t,] <- resample(w[t-1,])
      logZ = logZ + normalise_weights_in_log_space(w[t-1,])[[2]]
      
      for(i in 1:N){
        #filtering particles
        X[t,i,] <- sample_twisted_transition(X[t-1, ancestors[t,i],], model, psi_pa[t,], 1)
        w[t,i] <- eval_twisted_potential(model, list(NA, NA, psi_pa[t,]), X[t,i,], obs[t,])
        
      }
    }else{
      ancestors[t,] <- ancestors[t-1,]
      for(i in 1:N){
        
        X[t,i,] <- sample_twisted_transition(X[t-1, i,], model, psi_pa[t,], 1)
        w[t,i] <- w[t-1,i] + eval_twisted_potential(model, list(NA, NA, psi_pa[t,]), X[t,i,], obs[t,])
      }
    }
    
    logZ <- logZ + normalise_weights_in_log_space(w[Time,])[[2]] 
  }
  
  
  
  
  return(list(X, w, logZ, ancestors))
}
#' @import mvnfast
