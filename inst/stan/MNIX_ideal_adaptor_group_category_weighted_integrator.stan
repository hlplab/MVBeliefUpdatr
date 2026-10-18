/*
  MNIX Ideal Adaptor - Group-Category-Weighted Integrator

  Fit multinomial response data using a belief-updating model to infer prior
  parameters and group-by-category attention/cue weights under late cue integration.

  Cognitive & Biological Interpretation:
  Under late integration with a group-category-weighted integrator, attention weights
  w_{l, c} in Delta^{K-1} vary across both exposure groups l and categories c:
    log p(x | c) = sum_{k=1}^K w_{l, c, k} log p(x_k | c) + log prior_p[c]
  This full tensor model allows the perceptual weighting of cues to differ across
  categories while also adapting empirically across exposure conditions.

  Mathematical Note on Coordinate Scaling Invariance:
  Attention weights w_{l, c, k} act as precision allocations across marginal
  log-densities. Under diagonal affine scaling x_k = a_k * x_{k,orig} + b_k (a_k > 0),
  the additive constants reflect coordinate scalings. The dimensionless precision
  allocations remain invariant and do not require backtransformation.
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

  int<lower=0, upper=1> weights_known;                     // 1 if weights fixed by user, 0 if estimated
  int<lower=0, upper=1> lapse_rate_known;
  int<lower=0, upper=1> lapse_bias_known;
  int<lower=0, upper=1> mu_0_known;
  int<lower=0, upper=1> Sigma_0_known;
  array[weights_known ? L : 0, weights_known ? M : 0] simplex[K] weights_data; // optional: user-provided weights
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
  array[weights_known ? 0 : L, weights_known ? 0 : M] simplex[K] weights_param; // group-category attention weights
  real<lower=K> kappa_0;                                        // prior pseudocount for category mu
  real<lower=K+1> nu_0;                                         // prior pseudocount for category Sigma

  array[lapse_rate_known ? 0 : 1] real<lower=0, upper=1> lapse_rate_param;

  array[mu_0_known ? 0 : M] vector[K] m_0_param_transformed;   // prior expected means (standardized)
  vector<lower=0>[mu_0_known ? 0 : K] m_0_tau_transformed;     // prior SD of m_0
  array[Sigma_0_known ? 0 : M] vector<lower=0>[K] tau_0_param_transformed; // prior scale of Chisquare (standardized)
}

transformed parameters {
  array[L, M] simplex[K] weights = weights_known ? weights_data : weights_param;
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

  // compute posterior probabilities for each test stimulus (group-category weighted integration)
  for (j in 1:n_test) {
    int group = y_test[j];
    for (cat in 1:M) {
      vector[K] cue_log_lik;
      for (k in 1:K) {
        cue_log_lik[k] = student_t_lpdf(x_test[j][k] |
                                        nu_n[cat,group],
                                        m_n_transformed[cat,group][k],
                                        t_scale_transformed[cat,group][k]);
      }
      log_p_test_conj[j,cat] = dot_product(weights[group, cat], cue_log_lik) + log(prior_p[cat]);
    }
    p_test_conj[j] = exp(log_p_test_conj[j] - log_sum_exp(log_p_test_conj[j]));
  }
}

model {
  if (!weights_known) {
    for (g in 1:L) {
      for (cat in 1:M) {
        weights_param[g, cat] ~ dirichlet(rep_vector(1.0, K));
      }
    }
  }
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

  for (cat in 1:M) {
    m_0[cat] = unscale_mean(INV_SCALE, shift, m_0_transformed[cat]);
    S_0[cat] = unscale_diag_var(INV_SCALE, tau_0_transformed[cat]);
    for (group in 1:L) {
      m_n[cat,group] = unscale_mean(INV_SCALE, shift, m_n_transformed[cat,group]);
      S_n[cat,group] = unscale_diag_var(INV_SCALE, tau_n_transformed[cat,group]);
    }
  }
}
