#ess resampling

#logW = w[t-1,]

ESS <- function(logW, is.log=FALSE){
  if(is.log) {
    mx <- max(logW)
    s <- sum(exp(logW - mx))
    ess <- 1/sum((exp(logW - mx)/s)^2)
  }else{
    s <- sum(logW)
    ess <- 1/sum((logW/s)^2) 
  }
  return(ess)  
}