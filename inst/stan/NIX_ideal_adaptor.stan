/*
  NIX Ideal Adaptor (1D)

  Fit multinomial response data using a belief-updating model to infer prior
  parameters under 1D Normal-Inverse-Chi-Squared (NIX) adaptation.
*/

functions {
#include include/mvbu_functions.stan
}


data {
  int M;                              // number of categories
  int L;                              // number of groups (conditions or subjects)

  array[M,L] int<lower=0> N_exposure; // observations per category (M) and group (L)
  matrix[M,L] x_mean_exposure;        // mean of training observations
  matrix[M,L] x_css_exposure;         // sum of centered squares of training observations

  int N_test;                         // number of test trials
  array[N_test] real x_test;          // locations of test trials
  array[N_test] int y_test;           // group labels for test trials
  array[N_test,M] int z_test_counts;  // responses for test trials

  simplex[M] prior_p;                 // prior category probabilities

  int<lower=0, upper=1> lapse_rate_known;
  int<lower=0, upper=1> lapse_bias_known;
  int<lower=0, upper=1> mu_0_known;
  int<lower=0, upper=1> Sigma_0_known;
  array[lapse_rate_known ? 1 : 0] real<lower=0, upper=1> lapse_rate_data;
  array[lapse_bias_known ? 1 : 0] simplex[M] lapse_bias_data;
  array[mu_0_known ? M : 0] real mu_0_data;
  array[Sigma_0_known ? M : 0] real Sigma_0_data;

  real<lower=0> tau_scale;            // scale for the prior over means and variances

  real INV_SCALE;                     // Inverse of transformation applied to data
  real shift;                         // Center of affine transformation applied to data

  int<lower=0, upper=1> split_loglik_per_observation;
}

transformed data {
  simplex[M] lapse_bias = lapse_bias_known ? lapse_bias_data[1] : rep_vector(1.0 / M, M);
  real<lower=0> sigma_kappanu = max(max(to_array_1d(N_exposure)) * 4, 10);                // scale for the prior of kappa/nu_0 (at least 10)
  real m_0_mu = 0;                                                                        // center of prior of m_0
}

parameters {
  /* The strength of the prior beliefs is currently assumed to be shared across groups and categories.
     Any of these assumptions could be revisited in the future. For example, different groups could differ in
     the strength of their prior beliefs because of differences in prior experience. Different categories could
     be associated with weaker or stronger prior beliefs because of differences in the amount of exposure
     (overall frequency of the category) or because differences in e.g., the cross-talker variability in the
     realization of the category.
  */
  real<lower=1> kappa_0;             // prior pseudocount for mean
  real<lower=2> nu_0;                // prior pseudocount for sd

  array[lapse_rate_known ? 0 : 1] real<lower=0, upper=1> lapse_rate_param;
  array[mu_0_known ? 0 : M] real m_0_param_transformed;                  // prior expected means (standardized)
  array[mu_0_known ? 0 : 1] real<lower=0> m_0_tau_transformed;           // prior sd of expected means
  array[Sigma_0_known ? 0 : M] real<lower=0> tau_0_param_transformed;    // prior sd of scatter (standardized)
}

transformed parameters {
  array[M] real m_0_transformed = mu_0_known ? mu_0_data : m_0_param_transformed;
  array[M] real S_0_transformed;
  real lapse_rate = lapse_rate_known ? lapse_rate_data[1] : lapse_rate_param[1];

  // updated beliefs depend on input/group
  array[M,L] real<lower=1> kappa_n;                // updated mean pseudocount
  array[M,L] real<lower=2> nu_n;                   // updated sd pseudocount
  array[M,L] real m_n_transformed;                 // updated expected mean
  array[M,L] real<lower=0> S_n_transformed;        // updated expected scatter
  array[M,L] real<lower=0> t_scale_transformed;    // scale parameter of predictive t distribution

  array[N_test] simplex[M] p_test_conj;
  array[N_test] vector[M] log_p_test_conj;

  // update NIX2 parameters according to conjugate updating rules. Note that we're using the 
  // scatter parameterization (S_0, S_n) rather than the variance/scale parameterization of 
  // Murphy (2012). This keeps things parallel to the NIW model, and also means that post-
  // processing functions like get_expected_sigma_from_S have the same meaning for both 
  // types of models (rather than double transforming the NIX scale parameter when the goal 
  // is to extract the expected category variance).
  for (cat in 1:M) {
    S_0_transformed[cat] = Sigma_0_known ? Sigma_0_data[cat] * (nu_0 - 2) : tau_0_param_transformed[cat]^2;
    for (group in 1:L) {
      kappa_n[cat,group] = kappa_0 + N_exposure[cat,group];
      nu_n[cat,group] = nu_0 + N_exposure[cat,group];
      m_n_transformed[cat,group] =
        nix_1d_update_mean(kappa_0, m_0_transformed[cat], N_exposure[cat,group], x_mean_exposure[cat,group]);
      S_n_transformed[cat,group] =
        nix_1d_update_S(kappa_0, S_0_transformed[cat], N_exposure[cat,group], x_mean_exposure[cat,group], x_css_exposure[cat,group], m_0_transformed[cat]);
      t_scale_transformed[cat,group] =
        nix_1d_predictive_scale(S_n_transformed[cat,group], kappa_n[cat,group], nu_n[cat,group]);
    }
  }

  // compute category probabilities for each test stimulus
  for (j in 1:N_test) {
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
      m_0_param_transformed ~ normal(m_0_mu, m_0_tau_transformed[1]);
  }

  if (!Sigma_0_known) {
    for (cat in 1:M) {
      tau_0_param_transformed[cat] ~ cauchy(0, tau_scale);
    }
  }

  for (n in 1:N_test) {
    z_test_counts[n] ~ multinomial(lapse_mixture(p_test_conj[n], lapse_rate, lapse_bias));
  }
}

generated quantities {
#include include/gq_multinomial_loglik.stan

  // Backtransformed parameters in natural physical units (clean base names)
  array[M] real m_0;
  array[M] real S_0;
  array[M,L] real m_n;
  array[M,L] real S_n;
  for (cat in 1:M) {
    m_0[cat] = unscale_1d_mean(INV_SCALE, shift, m_0_transformed[cat]);
    S_0[cat] = unscale_1d_scatter(INV_SCALE, S_0_transformed[cat]);
    for (group in 1:L) {
      m_n[cat,group] = unscale_1d_mean(INV_SCALE, shift, m_n_transformed[cat,group]);
      S_n[cat,group] = unscale_1d_scatter(INV_SCALE, S_n_transformed[cat,group]);
    }
  }
}
