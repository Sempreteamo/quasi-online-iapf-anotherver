normalise_weights_in_log_space <- function(logW_un){
  mx <- max(logW_un)
  logW <- exp(logW_un - mx)/sum(exp(logW_un - mx))
  logZ <- log(mean(exp(w_apf[n,]-mx))) + mx
  return(list(logW, logZ))
}