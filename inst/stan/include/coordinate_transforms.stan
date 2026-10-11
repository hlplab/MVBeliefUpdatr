/*
  Coordinate Transformation and Unscaling Functions for MVBeliefUpdatr
  
  Covers:
  - Backtransforming standardized parameters to original physical cue coordinates
  - Vector, diagonal variance, full covariance matrix, and 1D transformations
*/

// Unscale mean vector: INV_SCALE * m_std - shift
vector unscale_mean(matrix INV_SCALE, vector shift, vector m_std) {
  return INV_SCALE * m_std - shift;
}

// Unscale diagonal category variance vector from standard deviations tau_std
vector unscale_diag_var(matrix INV_SCALE, vector tau_std) {
  return diagonal(INV_SCALE * diag_matrix(tau_std^2) * INV_SCALE');
}

// Unscale full category covariance matrix from scatter matrix S_std
matrix unscale_cov(matrix INV_SCALE, matrix S_std) {
  return INV_SCALE * S_std * INV_SCALE';
}

// Unscale 1D mean on composite/single feature: s_w * m_std - shift_w
real unscale_1d_mean(real s_w, real shift_w, real m_std) {
  return s_w * m_std - shift_w;
}

// Unscale 1D variance from scatter S_std and degrees of freedom nu
real unscale_1d_var(real s_w, real S_std, real nu) {
  return (s_w^2 * S_std) / (nu - 2.0);
}

// Unscale 1D scatter on composite/single feature: s_w^2 * S_std
real unscale_1d_scatter(real s_w, real S_std) {
  return s_w^2 * S_std;
}

