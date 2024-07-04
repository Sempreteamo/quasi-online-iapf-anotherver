compute_dKS <- function(x, w, fks.obj){
  d <-  dim(x)[2]
  w_ <- normalise_weights_in_log_space(w)
  weighted_mean <- colSums(w_*x)
  
  fks_mean <- fks.obj$ahatt[,Time]
  fks_cov <- fks.obj$Vt[,,Time]
  
  dist <- mahalanobis(weighted_mean, fks_mean, fks_cov)/d
  
  return(dist)
}
