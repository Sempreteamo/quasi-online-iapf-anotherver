# params <- list(dist_mu, dist_cov)
# psi <- psi_pa[t,]
compute_twisted_params <- function(params, psi){
  d <- length(psi)/2
  dist_mu <- params$mean
  dist_cov <- params$cov
  psi_mu <- psi[1:d]
  psi_cov <- psi[(d+1):(d+d)]
  
  zc <- (psi_cov^(-1) + diag(dist_cov)^(-1))^(-1)

  mu <- zc*(diag(dist_cov)^(-1)*dist_mu + psi_cov^(-1)*psi_mu)
  
  cov <- diag(zc, nrow = d, ncol = d)
  
  return(params = list(mu, cov))
}

#test: 
#psi <- rnorm(2*d)
#compute_twisted_initial_params(params, psi)
