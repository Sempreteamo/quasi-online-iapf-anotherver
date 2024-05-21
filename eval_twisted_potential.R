#if the time t at psi_pa[t,] is 1, then density should be 0
#if time t at psi_pa[t,] is T, set the psi_pa of psi_til as NA
#if t is T+1, set psi_pa of psi as NA

eval_twisted_potential <- function(model, psi_pa, x, y){
  ini_mu <- model$ini_mu
  ini_cov <- model$ini_cov
  psi_d <- psi_pa[[1]]
  psi_t <- psi_pa[[2]]
  psi <- psi_pa[[3]]
  
  dif <- ini_mu - psi_d[1:d]
  add <- psi_d[(d+1):(d+d)] + diag(ini_cov)
  
    density <- (d/2)*log(2*pi)-(1/2)*log(det(diag(add, nrow=d,ncol=d))) - 
      (1/2)*t(dif)%*%diag((add)^(-1), nrow=d,ncol=d)%*%(dif)
    
    if(is.na(density)){
      density <- 0
    }
    
    potential <- evaluate_log_g(model, x, y) + evaluate_psi_tilde(x, psi_t, model) - density -
      evaluate_psi(x, psi) #g_2:T 
  
  
  return(potential)
}

#test init_mu <- rnorm(d)
#psi_pa <- rnorm(2*d)
#psi_pa <- rep(NA, d)
#y <- rnorm(d)
