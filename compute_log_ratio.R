log_ratio <- function(Z, model, data){
  a0 <- model$a0
  P0 <- model$P0
  dt <- model$dt
  ct <- model$ct
  Tt <- model$Tt
  Zt <- model$Zt
  Ht <- model$Ht
  Gt <- model$Gt
  obs <- data$obs
  fkf.obj <- fkf(a0, P0, dt, ct, Tt, Zt, Ht, Gt, yt = t(obs))
  fks.obj <- fks(fkf.obj)
  fkf.obj_Z <- fkf(a0, P0, dt, ct, Tt, Zt, Ht, Gt, yt = t(obs))$logLik
  cat('NC = ', fkf.obj_Z)
  cat('log_ratio = ', exp(Z-fkf.obj_Z))
}