NULL

# Recommended by Ben Bolker for efficiency
# from https://stats.stackexchange.com/questions/62850/obtaining-covariance-matrix-from-correlation-matrix
cor2cov = function(omega, tau) {
  outer(tau,tau) * omega
}

cov2tau = function(v) {
  sqrt(diag(v))
}

