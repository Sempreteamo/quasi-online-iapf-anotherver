#params <- list(ini_mu, ini_cov)
#psi <- psi_pa[t,]
sample_twisted_initial <- function(model, psi, N){
  
  output <- compute_twisted_params(model, psi)
  
  mu <- output[[1]]
  
  cov <- output[[2]]
  #set.seed(1234)
  samples <- rmvn(N, mu, cov)
  
  return(samples)
}

#test
#psi <- rnorm(2*d)
