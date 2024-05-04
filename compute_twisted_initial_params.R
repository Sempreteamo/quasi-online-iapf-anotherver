# params <- list(ini_mu, ini_cov)
# psi <- psi_pa[t,]
compute_twisted_initial_params <- function(params, psi){
  ini_mu <- params[[1]]
  ini_cov <- params[[2]]
  psi_mu <- psi[1:d]
  psi_cov <- psi[(d+1):(d+d)]
  
  zc <- (psi_cov^(-1) + diag(ini_cov)^(-1))^(-1)
  
  mu <- zc*(diag(ini_cov)^(-1)*ini_mu + psi_cov^(-1)*psi_mu)
  
  cov <- diag(zc, nrow = d, ncol = d)
  
  return(twisted_params = list(mu, cov))
}

#test: psi <- rnorm(2*d)
