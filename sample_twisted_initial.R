sample_twisted_initial <- function(params, psi, N){
  
  
  output <- compute_twisted_initial_params(params, psi)
  
  mu <- output[[1]]
  
  cov <- output[[2]]
  
  samples <- rmvn(N, mu, cov)
  
  return(samples)
}
