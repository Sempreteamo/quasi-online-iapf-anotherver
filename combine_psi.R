combine_psi <- function(psi, index){
  Time <- nrow(psi[[1]])
  
  combined_psi <- matrix(NA, Time, 2*d)
    
   for(t in 1: Time){
     combined_psi[t,] <- psi[[index[t]]][t,]
   }
    
    return(psi = combined_psi)
}

#test: 
#psi1 <- matrix(1:10, 10, 1)
#psi2 <- matrix(11:20, 10, 1)
#psi <- list(psi1, psi2)
#index <- sample(1:2, 10, replace = TRUE)
#combine_psi(psi, index)
