combine_psi <- function(psi_pa, data){
  index <- data$psi_index
  obx <- data$obs
  Time <- nrow(obs)
  
  combined_psi <- matrix(NA, Time, 2*d)
    
   for(t in 1: Time){
     combined_psi[t,] <- psi_pa[[index[t]]][t,]
   }
    
    return(psi = combined_psi)
}

#test: psi_pa1 <- matrix(rnorm(Time*2*d), Time, 2*d)
#psi_pa2 <- matrix(rnorm(Time*2*d), Time, 2*d)
#psi_pa <- list(psi_pa1, psi_pa2)
