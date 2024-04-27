#logW = w[t-1,] or w[t-1,]
#Z_apf = Z_apf[l]
#N = N[l]
resample <- function(logW, mode = 'res'){
  mx <- max(logW)
  w_ <- exp(logW - mx)/sum(exp(logW - mx))
  
  if(mode == 'multi'){
    indices <- sample(1:N, N, replace = TRUE, prob = w_)
    return(indices)
    
  }else{
    Ntm <- as.integer(Num*w_)
    
    indices <- unlist(lapply(1:Num, function(i) {rep(i, Ntm[i])}))
    mr <- Num - sum(Ntm)
    
    w_hat <- w_ - Ntm/Num
    w_hat <- w_hat*Num/mr
    
    indices <- c(sample(1:Num, mr, replace = TRUE, prob = w_hat), indices)
    return(indices)
  }
}