
#params <- list(A, B)
#psi <- psi_pa[t,]

sample_twisted_transition <- function(x, params, psi, N){
  A <- params[[1]]%*%x
  
  params <- list(A, params[[2]])
  
  output <- compute_twisted_initial_params(params, psi)
  
  mu <- output[[1]]
  
  cov <- output[[2]]
  
  param <- list(mu, cov)
  
  samples <- sample_normal_distribution(param, N)
  
  return(samples)
}

#test: x <- rnorm(d)
#psi <- rnorm(2*d)
