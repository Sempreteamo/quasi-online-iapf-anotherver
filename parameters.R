Napf = N = 1000
lag = 10
Time = 200
d = 2
alpha = 0.42
A <- matrix(nrow = d, ncol = d)
for (i in 1:d){
  for (j in 1:d){
    A[i,j] = alpha^(abs(i-j) + 1)
  }
}
ini_mu <- rep(0, d)
B = C = D = ini_cov = diag(1, nrow = d, ncol = d)
k <- 5
tau <- 0.5
kappa = 0.5
model <- list(ini_mu, ini_cov, A, B, C, D, k, tau, kappa)


data <- generate_obs()


