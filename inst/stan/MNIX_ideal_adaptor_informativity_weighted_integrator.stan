/*
  MNIX Ideal Adaptor - Informativity-Weighted Integrator

  Fit multinomial response data using a belief-updating model with dynamic
  informativity weighting based on Toscano & McMurray (2010).

  Cognitive & Biological Interpretation:
  In normative Bayesian cue combination (Ernst & Banks 2002; Landy et al. 2011),
  an ideal observer uses Bayes rule, where the relative informativity of cues
  is already fully reflected in the likelihood ratio defined by the adapted category
  distributions (ideal integrator).
  In contrast, Toscano & McMurray (2010) ("Cue weighting with competition: A neural
  network model of speech perception", Cognitive Science) proposed an empirical
  heuristic model where perceptual attention weights w_{l, k} adapt dynamically
  based on between-category cue discriminability:
    I_{l, k} = sum_{c1 < c2} prior_p[c1] * prior_p[c2] * (m_{l, c1, k} - m_{l, c2, k})^2 / (tau_{l, c1, k} * tau_{l, c2, k})
    w_{l} = softmax(I_{l})
  Because this model applies a multiplicative attention weight w_{l, k} on top of
  the adapted likelihoods, it "double-counts" cue informativity relative to Bayes rule.
  It represents a dynamic informativity-weighted heuristic model, NOT a normative ideal integrator.

  Mathematical Note on Coordinate Scaling Invariance:
  Dynamic weights w_{l, k} are computed as a softmax over standardized discriminabilities,
  acting as dimensionless precision allocations across marginal Student-t log-densities.
  The denominator factorization exp(-C) ensures that the resulting posterior categorization
  probabilities are invariant to coordinate scaling and translation.
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

  array[L] simplex[K] cue_weight_n; // dynamically adapted informativity weights per group

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

  // Compute Toscano & McMurray (2010) informativity weights dynamically for each group
  for (group in 1:L) {
    vector[K] raw_weight = rep_vector(0, K);
    for (cat1 in 1:M) {
      for (cat2 in 1:M) {
        raw_weight += prior_p[cat1] * prior_p[cat2] *
          ((m_n_transformed[cat1,group] - m_n_transformed[cat2,group])^2 ./
           (tau_n_transformed[cat1,group] .* tau_n_transformed[cat2,group]));
      }
    }
    raw_weight ./= rep_vector(2.0, K);
    cue_weight_n[group] = softmax(raw_weight);
  }

  // compute posterior probabilities for each test stimulus (informativity-weighted integration)
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
      log_p_test_conj[j,cat] = dot_product(cue_weight_n[group], cue_log_lik) + log(prior_p[cat]);
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
  array[L] simplex[K] weights = cue_weight_n;

  for (cat in 1:M) {
    m_0[cat] = unscale_mean(INV_SCALE, shift, m_0_transformed[cat]);
    S_0[cat] = unscale_diag_var(INV_SCALE, tau_0_transformed[cat]);
    for (group in 1:L) {
      m_n[cat,group] = unscale_mean(INV_SCALE, shift, m_n_transformed[cat,group]);
      S_n[cat,group] = unscale_diag_var(INV_SCALE, tau_n_transformed[cat,group]);
    }
  }
}
