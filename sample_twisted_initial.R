#params <- list(ini_mu, ini_cov)
#psi <- psi_pa[t,]
sample_twisted_initial <- function(params, psi, N){
  
  output <- compute_twisted_params(params, psi)
  
  mu <- output[[1]]
  
  cov <- output[[2]]
  
  samples <- rmvn(N, mu, cov)
  
  return(samples)
}

#test
#psi <- rnorm(2*d)
