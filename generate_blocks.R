generate_blocks <- function(lag, len){
  num_blocks <- ceiling(len / lag)
  
  psi_index <- numeric(len)
  
  breaks <- list()
  
  breaks[[1]] <- seq(1, len, by = lag)
  breaks[[2]] <- c(1, seq(lag/2 + 1, len, by = lag))
  
  for (i in 1:num_blocks) {
    
    start_time <- (i - 1) * lag + 1
    end_time <- min(i * lag, len)
    
    if (i == 1) {
      start_time <- start_time
      end_time <- floor((3/4) * lag)
    }
    else {
      start_time <- start_time + floor((1/4) * lag)
      end_time <- end_time - floor((1/4) * lag)
    }
    
    psi_index[start_time:end_time] <- 1
  }
  
  psi_index[psi_index == 0] <- 2
  
  return(list(breaks, psi_index))
}