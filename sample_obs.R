#param <- list(A, B, C, D)
#N = Time
sample_obs <- function(param, N){
  A <- param[[1]]
  B <- param[[2]]
  C <- param[[3]]
  D <- param[[4]]
  
  X <-  matrix(0, nrow = N, ncol = d)
  
  X[1,] <- rnorm(d)     
  
  for(t in 2:Time){ 
    X[t,] <- chol(B)%*%rnorm(d) + A%*%X[t-1,]  #t(rmvn(d) + A%*%x)
  }
  data <- sample_g(list(t(C%*%t(X)), D), X, N)
  return(data)
}
