#logW = w[t-1,] or w[t-1,]
#Z_apf = Z_apf[l]
#N = N[l]
resample <- function(logW, mode = 'res'){
  
  N <- length(logW)
  w_ <- normalise_weights_in_log_space(logW)[[1]]
  
  if(mode == 'multi'){
    indices <- sample(1:N, N, replace = TRUE, prob = w_)
    return(indices)
    
  }else{
    Ntm <- as.integer(N*w_)
    
    indices <- unlist(lapply(1:N, function(i) {rep(i, Ntm[i])}))
    mr <- N - sum(Ntm)
    
    w_hat <- w_ - Ntm/N
    w_hat <- w_hat*N/mr
    
    indices <- c(sample(1:N, mr, replace = TRUE, prob = w_hat), indices)
    return(indices)
  }
}
