evaluate_log_g <- function(y, x){  
  log_density <- (-d/2)*log(2*pi) - (1/2)*log(det(D)) + (-1/2)*t(y-C%*%x)%*%D%*%(y-C%*%x)
  return (log_density) #obs prob  C%*%x = x
}