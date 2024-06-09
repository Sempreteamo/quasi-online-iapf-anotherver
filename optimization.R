#x = X_apf[t,,]
#lfn = log(psi[t,])

optimization <- function(x, lfn){
  params <- vector()
  d <- dim(x)[2]
  
  coef <- -lm(lfn~., data.frame(cbind(x^2, x)))$coefficients
  a <- coef[2:(1+d)]
  b <- coef[(2+d):length(coef)]
  
  params <- c(b/(-2*a), 1/(2*a))
  
  return(params)
}
