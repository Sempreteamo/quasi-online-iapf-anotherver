#param <- list(A, B, C, D)
#N = Time
sample_obs <- function(model, N){
  A <- model$A
  B <- model$B
  C <- model$C
  D <- model$D
  
  X <- matrix(0, nrow = N, ncol = d)
  data <- matrix(0, nrow = N, ncol = d)
  
  X[1,] <- rnorm(d)     
  
  for(t in 2:Time){ 
    X[t,] <- chol(B)%*%rnorm(d) + A%*%X[t-1,]  #t(rmvn(d) + A%*%x)
  }
  for(t in 1:N){
    data[t,] <- rmvn(1, C%*%X[t,], D)
  }
  return(data)
}
