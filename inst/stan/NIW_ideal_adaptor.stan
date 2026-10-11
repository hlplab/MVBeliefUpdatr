/*
* Fit multinomial response data using a belief-updating model to infer prior
* parameters.  A normal-inverse-Wishart prior is used, and it's assumed
* that the subject knows the true labels of all the input stimuli (e.g.,
* doesn't model any uncertainty in classification).
*
* Sufficient statistics and (if provided) prior mu_0 and Sigma_0 are assumed
* to be in the space defined by an affine transformation f(x) = SCALE(x + shift).
* Typically, that's the space obtained by whitening or standardizing the exposure
* data.
*
* Modified by
* Florian Jaeger
* April 2025
*/

functions {
#include include/mvbu_functions.stan
}


data {
  int K;                                                   // number of features
  int L;                                                   // number of exposure groups (e.g. subjects)
  int M;                                                   // number of categories

  array[M,L] int<lower=0> N_exposure;                      // number of observations per category (M) and exposure group (L)
  array[M,L] vector[K] x_mean_exposure;                    // means for each category (M) and exposure group (L)
  array[M,L] matrix[K,K] x_css_exposure;                   // sum of centered squares matrix for each category (M) and group (L)
                                                           // (using matrix, rather than cov_matrix to allow non-positive-definiteness, 
                                                           // necessary for groups without any exposure, which have matrices of 0s)

  int N_test;                                              // number of unique combinations of test locations & exposure groups
  array[N_test] vector[K] x_test;                          // locations (in cue space) of test trials
  array[N_test] int y_test;                                // exposure group indicator of test trials
  array[N_test,M] int z_test_counts;                       // responses for test trials (for each category M)

  simplex[M] prior_p;                                      // prior category probabilities

  int<lower=0, upper=1> lapse_rate_known;
  int<lower=0, upper=1> lapse_bias_known;
  int<lower=0, upper=1> mu_0_known;
  int<lower=0, upper=1> Sigma_0_known;
  array[lapse_rate_known ? 1 : 0] real<lower=0, upper=1> lapse_rate_data;       // optional: user provided lapse_rate
  array[lapse_bias_known ? 1 : 0] simplex[M] lapse_bias_data;                   // optional: user provided lapse bias
  array[mu_0_known ? M : 0] vector[mu_0_known ? K : 0] mu_0_data;               // optional: user provided expected mu_0 (prior category means) in space of affine transformation defined by INV_SCALE^-1 and shift
  array[Sigma_0_known ? M : 0] cov_matrix[Sigma_0_known ? K : 0] Sigma_0_data;  // optional: user provided expected Sigma_0 (prior category covariance matrices) in space of affine transformation defined by INV_SCALE^-1 and shift

  /* For now, this script assumes that the observations (cue vectors) are centered. The prior
     mean of m_0 is set to 0. In the future, one could automatically adjust the location based
     on the input data (e.g., by placing the mean of the Normal prior in the center of the
     exposure data).
  */
  vector<lower=0>[K] tau_scale;                            // separate taus for each of the K features to capture that features can be on separate scales
  real<lower=0> L_omega_eta;                               // scale of LKJ prior for correlation of variance of m_0

  /* The data above are assumed to have been transformed by a sensible affine transformation f(x) = A(x - center)
     (e.g., by scaling or whitening the data) in order to improve numerical stability. A_inv and center are used
     below in order to transform the parameters of interest back into the original space.
  */
  matrix[K, K] INV_SCALE;                                  // Inverse of transformation matrix (must be invertible) that was applied to data
  vector[K] shift;                                         // Center of affine transformation that was applied to data

  int<lower=0, upper=1> split_loglik_per_observation;
}

transformed data {
  simplex[M] lapse_bias = lapse_bias_known ? lapse_bias_data[1] : rep_vector(1.0 / M, M);
  real<lower=0> sigma_kappanu = max(max(to_array_1d(N_exposure)) * 4, 10); // scale for the prior of kappa/nu_0 (at least 10)
  vector[K] m_0_mu = rep_vector(0, K);                     // center of prior of m_0
}

parameters {
  /* The strength of the prior beliefs is currently assumed to be shared across groups, categories, and cues.
     Any of these assumptions could be revisited in the future. For example, different groups could differ in
     the strength of their prior beliefs because of differences in prior experience. Different categories could
     be associated with weaker or stronger prior beliefs because of differences in the amount of exposure
     (overall frequency of the category) or because differences in e.g., the cross-talker variability in the
     realization of the category. Similar reasoning can be applied to cues.
  */
  real<lower=K> kappa_0;                                   // prior pseudocount for category mu
  real<lower=K+1> nu_0;                                    // prior pseudocount for category Sigma

  array[lapse_rate_known ? 0 : 1] real<lower=0, upper=1> lapse_rate_param;

  array[mu_0_known ? 0 : M] vector[K] m_0_param_transformed;           // prior expected means (standardized)
  vector<lower=0>[mu_0_known ? 0 : K] m_0_tau_transformed;             // prior variances of expected means (m_0)
  cholesky_factor_corr[mu_0_known ? 0 : K] m_0_L_omega_transformed;    // prior correlations of variances of expected means (in cholesky form)

  array[Sigma_0_known ? 0 : M] vector<lower=0>[K] tau_0_param_transformed;          // standard deviations of prior scatter matrix S_0
  array[Sigma_0_known ? 0 : M] cholesky_factor_corr[K] L_omega_0_param_transformed; // correlation matrix of prior scatter matrix S_0 (in cholesky form)
}

transformed parameters {
  array[M] vector[K] m_0_transformed = mu_0_known ? mu_0_data : m_0_param_transformed;                    // prior expected means m_0
  array[M] cov_matrix[K] S_0_transformed;                                                                  // prior scatter matrix S_0
  real lapse_rate = lapse_rate_known ? lapse_rate_data[1] : lapse_rate_param[1];                          // lapse rate

  // updated beliefs depend on input and group
  array[M,L] real<lower=K> kappa_n;          // updated mean pseudocount
  array[M,L] real<lower=K> nu_n;             // updated sd pseudocount
  array[M,L] vector[K] m_n_transformed;      // updated m_0
  array[M,L] cov_matrix[K] S_n_transformed;  // updated scatter matrix
  array[M,L] cov_matrix[K] t_scale_transformed; // scale matrix of predictive t distribution

  array[N_test] simplex[M] p_test_conj;
  array[N_test] vector[M] log_p_test_conj;

  /* Update NIW parameters according to conjugate updating rules taken from
     Murphy (2007, p. 134)
  */
  for (cat in 1:M) {
    // Get S_0 from expected Sigma given nu_0
    S_0_transformed[cat] = Sigma_0_known ? Sigma_0_data[cat] * (nu_0 - K - 1) : quad_form_diag(multiply_lower_tri_self_transpose(L_omega_0_param_transformed[cat]), tau_0_param_transformed[cat]);

    for (group in 1:L) {
      kappa_n[cat,group] = kappa_0 + N_exposure[cat,group];
      nu_n[cat,group] = nu_0 + N_exposure[cat,group];
      m_n_transformed[cat,group] =
        nix_update_mean(kappa_0, m_0_transformed[cat], N_exposure[cat,group], x_mean_exposure[cat,group]);
      S_n_transformed[cat,group] =
        niw_update_S(kappa_0, S_0_transformed[cat], N_exposure[cat,group], x_mean_exposure[cat,group], x_css_exposure[cat,group], m_0_transformed[cat]);
      t_scale_transformed[cat,group] =
        niw_predictive_scale(S_n_transformed[cat,group], kappa_n[cat,group], nu_n[cat,group], K);
    }
  }

  // compute posterior probabilities of all categories for each of the test stimuli
  for (j in 1:N_test) {
    int group = y_test[j];
    for (cat in 1:M) {
      // log likelihood of the test stimulus given the category
      log_p_test_conj[j,cat] = multi_student_t_lpdf(x_test[j] |
                                                    nu_n[cat,group] - K + 1,
                                                    m_n_transformed[cat,group],
                                                    t_scale_transformed[cat,group]) + log(prior_p[cat]);
    }
    // normalize to get actual posterior probs in simplex
    p_test_conj[j] = exp(log_p_test_conj[j] - log_sum_exp(log_p_test_conj[j]));
  }
}

model {
  kappa_0 ~ normal(0, sigma_kappanu);
  nu_0 ~ normal(0, sigma_kappanu);

  if (!mu_0_known) {
      m_0_tau_transformed ~ cauchy(0, tau_scale);
      m_0_L_omega_transformed ~ lkj_corr_cholesky(L_omega_eta);
      m_0_param_transformed ~ multi_normal_cholesky(m_0_mu, diag_pre_multiply(m_0_tau_transformed, m_0_L_omega_transformed));
  }

  // Specifying prior for components of S_0:
  if (!Sigma_0_known) {
    for (cat in 1:M) {
      tau_0_param_transformed[cat] ~ cauchy(0, tau_scale);
      L_omega_0_param_transformed[cat] ~ lkj_corr_cholesky(L_omega_eta);
    }
  }

  for (i in 1:N_test) {
    z_test_counts[i] ~ multinomial(lapse_mixture(p_test_conj[i], lapse_rate, lapse_bias));
  }

}

generated quantities {
#include include/gq_multinomial_loglik.stan

  // Backtransformed parameters in natural physical units (clean base names)
  array[M] vector[K] m_0;
  array[M] cov_matrix[K] S_0;
  array[M,L] vector[K] m_n;
  array[M,L] cov_matrix[K] S_n;
  for (cat in 1:M) {
    m_0[cat] = unscale_mean(INV_SCALE, shift, m_0_transformed[cat]);
    S_0[cat] = unscale_cov(INV_SCALE, S_0_transformed[cat]);
    for (group in 1:L) {
      m_n[cat,group] = unscale_mean(INV_SCALE, shift, m_n_transformed[cat,group]);
      S_n[cat,group] = unscale_cov(INV_SCALE, S_n_transformed[cat,group]);
    }
  }
}

