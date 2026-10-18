/*
  MNIX Ideal Adaptor - Ideal Integrator

  Fit multinomial response data using a belief-updating model to infer prior
  parameters under normative late cue integration.
  
  Cognitive & Biological Interpretation:
  Under late integration with an ideal integrator, cues are integrated according
  to unweighted conditional independence (pure Bayes rule):
    log p(x | c) = sum_{k=1}^K log p(x_k | c) + log prior_p[c]
  This serves as the normative Bayesian baseline. Changes in cue diagnosticity
  across exposure groups emerge adaptively and purely through belief updating of
  the marginal category distributions, with no free attention or weighting parameters.

  Mathematical Note on Coordinate Scaling Invariance:
  Under diagonal affine scaling of cues x_k = a_k * x_{k,orig} + b_k (a_k > 0):
  log p(x_k | c) = log p_orig(x_{k,orig} | c) - log |a_k|
  Summing across cues introduces an additive constant C = sum_k log |a_k| that has
  zero dependence on category c. In the posterior softmax probability P(c | x),
  the factor exp(-C) appears in both the numerator and every term of the denominator
  sum over all M categories. Factoring exp(-C) out of the denominator summation cancels
  it identically with the numerator. Thus, late integration is coordinate-invariant.
*/

functions {
  #include include/mvbu_functions.stan
}

data {
  int K;                                                   // number of features (cues >= 2)
  int L;                                                   // number of exposure groups (e.g. subjects)
  int M;                                                   // number of categories

  array[M,L] int<lower=0> n_exposure;                      // observations per category (M) and group (L)
  array[M,L] vector[K] mean_exposure;                    // means per category (M) and group (L)
  array[M,L] vector[K] css_exposure;                     // centered sum of squares per category and group

  int n_test;                                              // unique test locations & exposure groups
  array[n_test] vector[K] x_test;                          // stimulus locations in cue space
  array[n_test] int y_test;                                // group indicator of test trials
  array[n_test,M] int z_test_counts;                       // response counts per category M

  simplex[M] prior_p;                                      // prior category probabilities

  int<lower=0, upper=1> lapse_rate_known;
  int<lower=0, upper=1> lapse_bias_known;
  int<lower=0, upper=1> mu_0_known;
  int<lower=0, upper=1> Sigma_0_known;
  array[lapse_rate_known ? 1 : 0] real<lower=0, upper=1> lapse_rate_data;
  array[lapse_bias_known ? 1 : 0] simplex[M] lapse_bias_data;
  array[mu_0_known ? M : 0] vector[mu_0_known ? K : 0] mu_0_data;
  array[Sigma_0_known ? M : 0] cov_matrix[Sigma_0_known ? K : 0] Sigma_0_data;

  vector<lower=0>[K] tau_scale;                            // prior scale for cue SDs

  matrix[K, K] INV_SCALE;                                  // inverse affine transform matrix
  vector[K] shift;                                         // affine transform center

  int<lower=0, upper=1> split_loglik_per_observation;
}

transformed data {
  array[Sigma_0_known ? M : 0] vector[Sigma_0_known ? K : 0] tau_0_data;
  if (Sigma_0_known) {
    for (cat in 1:M) {
      tau_0_data[cat] = sqrt(diagonal(Sigma_0_data[cat]));
    }
  }
  simplex[M] lapse_bias = lapse_bias_known ? lapse_bias_data[1] : rep_vector(1.0 / M, M);
  real<lower=0> sigma_kappanu = max(max(to_array_1d(n_exposure)) * 4, 10);
  vector[K] m_0_mu = rep_vector(0, K);
}

parameters {
  real<lower=K> kappa_0;                                        // prior pseudocount for category mu
  real<lower=K+1> nu_0;                                         // prior pseudocount for category Sigma

  array[lapse_rate_known ? 0 : 1] real<lower=0, upper=1> lapse_rate_param;

  array[mu_0_known ? 0 : M] vector[K] m_0_param_transformed;   // prior expected means (standardized)
  vector<lower=0>[mu_0_known ? 0 : K] m_0_tau_transformed;     // prior SD of m_0
  array[Sigma_0_known ? 0 : M] vector<lower=0>[K] tau_0_param_transformed; // prior scale of Chisquare (standardized)
}

transformed parameters {
  array[M] vector[K] m_0_transformed = mu_0_known ? mu_0_data : m_0_param_transformed;
  array[M] vector<lower=0>[K] tau_0_transformed;
  real lapse_rate = lapse_rate_known ? lapse_rate_data[1] : lapse_rate_param[1];

  // updated beliefs depend on input and group
  array[M,L] real<lower=K> kappa_n;
  array[M,L] real<lower=K> nu_n;
  array[M,L] vector[K] m_n_transformed;
  array[M,L] vector[K] tau_n_transformed;
  array[M,L] vector[K] t_scale_transformed;

  array[n_test] simplex[M] p_test_conj;
  array[n_test] vector[M] log_p_test_conj;

  for (cat in 1:M) {
    tau_0_transformed[cat] = Sigma_0_known ? tau_0_data[cat] * sqrt((nu_0 - 2) / nu_0) : tau_0_param_transformed[cat];

    for (group in 1:L) {
      kappa_n[cat,group] = kappa_0 + n_exposure[cat,group];
      nu_n[cat,group] = nu_0 + n_exposure[cat,group];
      m_n_transformed[cat,group] = nix_update_mean(kappa_0, m_0_transformed[cat], n_exposure[cat,group], mean_exposure[cat,group]);
      tau_n_transformed[cat,group] = nix_update_tau(kappa_0, nu_0, tau_0_transformed[cat], n_exposure[cat,group], mean_exposure[cat,group], css_exposure[cat,group], m_0_transformed[cat]);
      t_scale_transformed[cat,group] = nix_predictive_scale(tau_n_transformed[cat,group], kappa_n[cat,group]);
    }
  }

  // compute posterior probabilities for each test stimulus (unweighted Bayes rule)
  for (j in 1:n_test) {
    int group = y_test[j];
    for (cat in 1:M) {
      log_p_test_conj[j,cat] = student_t_lpdf(x_test[j] |
                                              nu_n[cat,group],
                                              m_n_transformed[cat,group],
                                              t_scale_transformed[cat,group]) + log(prior_p[cat]);
    }
    p_test_conj[j] = exp(log_p_test_conj[j] - log_sum_exp(log_p_test_conj[j]));
  }
}

model {
  kappa_0 ~ normal(0, sigma_kappanu);
  nu_0 ~ normal(0, sigma_kappanu);

  if (!mu_0_known) {
    m_0_tau_transformed ~ cauchy(0, tau_scale);
    for (cat in 1:M) {
      m_0_param_transformed[cat] ~ normal(m_0_mu[cat], m_0_tau_transformed);
    }
  }

  if (!Sigma_0_known) {
    for (cat in 1:M) {
      tau_0_param_transformed[cat] ~ cauchy(0, tau_scale);
    }
  }

  for (i in 1:n_test) {
    z_test_counts[i] ~ multinomial(lapse_mixture(p_test_conj[i], lapse_rate, lapse_bias));
  }
}

generated quantities {
#include include/gq_multinomial_loglik.stan

  // Backtransformed parameters in natural physical units (clean base names)
  array[M] vector[K] m_0;
  array[M] vector[K] S_0;
  array[M,L] vector[K] m_n;
  array[M,L] vector[K] S_n;
  simplex[K] weights = rep_vector(1.0 / K, K);

  for (cat in 1:M) {
    m_0[cat] = unscale_mean(INV_SCALE, shift, m_0_transformed[cat]);
    S_0[cat] = unscale_diag_var(INV_SCALE, tau_0_transformed[cat]);
    for (group in 1:L) {
      m_n[cat,group] = unscale_mean(INV_SCALE, shift, m_n_transformed[cat,group]);
      S_n[cat,group] = unscale_diag_var(INV_SCALE, tau_n_transformed[cat,group]);
    }
  }
}
