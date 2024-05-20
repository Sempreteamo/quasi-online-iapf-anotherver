combine_psi <- function(psi_pa, data){
  index <- data$psi_index
  obs <- data$obs
  Time <- nrow(obs)
  
  combined_psi <- matrix(NA, Time, 2*d)
    
   for(t in 1: Time){
     combined_psi[t,] <- psi_pa[[index[t]]][t,]
   }
    
    return(psi = combined_psi)
}

#test: 
#psi_pa1 <- matrix(1:10, 10, 1)
#psi_pa2 <- matrix(11:20, 10, 1)
#psi_pa <- list(psi_pa1, psi_pa2)
#obs_ <- matrix(rnorm(10), 10, 1)
#index <- sample(1:2, 10, replace = TRUE)
#combine_psi(psi_pa, list(obs = obs_, psi_index = index))
