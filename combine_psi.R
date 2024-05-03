combine_psi <- function(psi, index){
    combined_psi <- matrix(NA, Time, d)
    
   for(t in 1: Time){
     combined_psi[t] <- psi[[index[t]]][t]
   }
    
    return(psi = combined_psi)
}