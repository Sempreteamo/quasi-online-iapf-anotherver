normalise_weights_in_log_space <- function(logW_un){
  mx <- max(logW_un)
  w_centered <- exp(logW_un - mx)
  sum_wc <- sum(w_centered)
  logW <- w_centered/sum_wc
  logZ <- log(sum_wc / length(logW_un)) + mx
  return(list(logW, logZ))
}
