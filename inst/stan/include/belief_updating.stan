/*
  Analytical Conjugate Belief Updating Functions for MVBeliefUpdatr
  
  Covers:
  - Vectorized NIX updates across K independent cues (late integration)
  - 1D scalar NIX updates on composite/feature dimensions
    (early integration and NIX)
  - Multivariate NIW updates for full covariance matrices
*/

// Vectorized NIX: updated mean
vector nix_update_mean(real kappa_0, vector m_0, int N, vector x_bar) {
  if (N > 0) {
    return (kappa_0 * m_0 + N * x_bar) / (kappa_0 + N);
  } else {
    return m_0;
  }
}

// Vectorized NIX: updated scale tau
vector nix_update_tau(
  real kappa_0, real nu_0, vector tau_0, int N,
  vector x_mean, vector x_css, vector m_0
) {
  if (N > 0) {
    real kappa_n = kappa_0 + N;
    real nu_n = nu_0 + N;
    vector[num_elements(tau_0)] ss_diff = (m_0 - x_mean)^2;
    return sqrt(
      (nu_0 * tau_0^2 + x_css + (N * kappa_0 / kappa_n) * ss_diff) / nu_n
    );
  } else {
    return tau_0;
  }
}

// Vectorized NIX: predictive Student-t scale
vector nix_predictive_scale(vector tau_n, real kappa_n) {
  return tau_n * sqrt((kappa_n + 1.0) / kappa_n);
}

// Vectorized NIX with known mean: updated scale tau
vector nix_known_mu_update_tau(
  real nu_0, vector tau_0, int N, vector x_mean, vector x_css, vector mu_0
) {
  if (N > 0) {
    real nu_n = nu_0 + N;
    vector[num_elements(tau_0)] ss_diff = (mu_0 - x_mean)^2;
    return sqrt((nu_0 * tau_0^2 + x_css + N * ss_diff) / nu_n);
  } else {
    return tau_0;
  }
}

// Vectorized NIX with known mean: predictive Student-t scale
vector nix_known_mu_predictive_scale(vector tau_n) {
  return tau_n;
}

// Vectorized NIX with known variance: predictive normal scale
vector nix_known_sigma_predictive_scale(vector tau_0, real kappa_n) {
  return tau_0 * sqrt((kappa_n + 1.0) / kappa_n);
}

// 1D Scalar NIX: updated mean
real nix_1d_update_mean(real kappa_0, real m_0, int N, real x_bar) {
  if (N > 0) {
    return (kappa_0 * m_0 + N * x_bar) / (kappa_0 + N);
  } else {
    return m_0;
  }
}

// 1D Scalar NIX: updated scatter S
real nix_1d_update_S(
  real kappa_0, real S_0, int N, real x_mean, real x_css, real m_0
) {
  if (N > 0) {
    real kappa_n = kappa_0 + N;
    return S_0 + x_css + ((N * kappa_0) / kappa_n) * (m_0 - x_mean)^2;
  } else {
    return S_0;
  }
}

// 1D Scalar NIX: predictive Student-t scale
real nix_1d_predictive_scale(real S_n, real kappa_n, real nu_n) {
  return sqrt((S_n * (kappa_n + 1.0)) / (kappa_n * nu_n));
}

// Multivariate NIW: updated scatter matrix S
matrix niw_update_S(
  real kappa_0, matrix S_0, int N, vector x_mean, matrix x_css, vector m_0
) {
  if (N > 0) {
    real kappa_n = kappa_0 + N;
    vector[num_elements(m_0)] diff_m = x_mean - m_0;
    return S_0 + x_css + ((N * kappa_0) / kappa_n) * (diff_m * diff_m');
  } else {
    return S_0;
  }
}

// Multivariate NIW: predictive Student-t scale matrix
matrix niw_predictive_scale(matrix S_n, real kappa_n, real nu_n, int K) {
  return S_n * (kappa_n + 1.0) / (kappa_n * (nu_n - K + 1.0));
}

// Multivariate NIW with known mean: updated scatter matrix S
matrix niw_known_mu_update_S(
  matrix S_0, int N, vector x_mean, matrix x_css, vector mu_0
) {
  if (N > 0) {
    vector[num_elements(mu_0)] diff_mu = x_mean - mu_0;
    return S_0 + x_css + N * (diff_mu * diff_mu');
  } else {
    return S_0;
  }
}

// Multivariate NIW with known mean: predictive Student-t scale matrix
matrix niw_known_mu_predictive_scale(matrix S_n, real nu_n, int K) {
  return S_n / (nu_n - K + 1.0);
}

// Multivariate normal with known covariance: predictive covariance matrix
matrix niw_known_sigma_predictive_scale(matrix Sigma_0, real kappa_n) {
  return Sigma_0 * (kappa_n + 1.0) / kappa_n;
}

// 1D Scalar NIX with known mean: updated scatter S
real nix_1d_known_mu_update_S(
  real S_0, int N, real x_mean, real x_css, real mu_0
) {
  if (N > 0) {
    return S_0 + x_css + N * (x_mean - mu_0)^2;
  } else {
    return S_0;
  }
}

// 1D Scalar NIX with known mean: predictive Student-t scale
real nix_1d_known_mu_predictive_scale(real S_n, real nu_n) {
  return sqrt(S_n / nu_n);
}

// 1D Scalar normal with known variance: predictive standard deviation
real nix_1d_known_sigma_predictive_scale(real sigma_0, real kappa_n) {
  return sigma_0 * sqrt((kappa_n + 1.0) / kappa_n);
}

