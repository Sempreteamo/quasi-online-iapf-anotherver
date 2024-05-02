#ess resampling

#logW_un = w[t-1,]

ESS_log <- function(logW_un){
    logW <- normalise_weights_in_log_space(logW_un)
    ess <- 1/sum(logW^2)
  return(ess)  
}
