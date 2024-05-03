
d = 3
alpha = 0.42
A <- matrix(nrow = d, ncol = d)
for (i in 1:d){
  for (j in 1:d){
    A[i,j] = alpha^(abs(i-j) + 1)
  }
}
ini_mu <- rep(0, d)
B = C = D = ini_cov = diag(1, nrow = d, ncol = d)
model <- list(ini_mu, ini_cov, A, B, C, D)


L = 10
Time = 200
lag <- generate_blocks(L, Time)[[1]]

data <- generate_obs()

Napf = N = 1000
