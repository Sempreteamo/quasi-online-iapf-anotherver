
#params <- list(A, B)
#psi <- psi_pa[t,]

sample_twisted_transition <- function(x, model, psi, N){
  A <- model$A
  B <- model$B
  
  params <- list(A%*%x, B)
  
  output <- compute_twisted_params(params, psi)
  
  mu <- output[[1]]
  
  cov <- output[[2]]
  
  samples <- rmvn(N, mu, cov)
  
  return(samples)
}
