/*
  MNIX Ideal Adaptor - Group-Weighted Integrator

  Fit multinomial response data using a belief-updating model to infer prior
  parameters and exposure group-specific attention/cue weights under late cue integration.

  Cognitive & Biological Interpretation:
  Under late integration with a group-weighted integrator, attention weights w_{l}
  in Delta^{K-1} vary across exposure groups l:
    log p(x | c) = sum_{k=1}^K w_{l, k} log p(x_k | c) + log prior_p[c]
  This captures empirical cue reweighting driven by experimental exposure conditions,
  such as shifts in perceptual attention induced by different training distributions.
  It provides a flexible empirical benchmark to test whether cue weighting changes
  differentially across exposure groups beyond normative belief updating.

  Mathematical Note on Coordinate Scaling Invariance:
  For each exposure group l, coordinate scaling x_k = a_k * x_{k,orig} + b_k adds
  C_l = sum_k w_{l, k} log |a_k|. In the posterior softmax probability P(c | x),
  the factor exp(-C_l) is constant across all categories within group l, factors
  outside the summation over categories in the denominator, and cancels identically.
  Thus, late integration group weights are invariant to coordinate scaling.
*/

functions {
  #include include/mvbu_functions.stan
}

data {
  int K; // number of features (cues >= 2)
  int L; // number of exposure groups (e.g. subjects)
  int M; // number of categories
  array[M, L] int<lower=0> n_exposure; // observations per category & group
  array[M, L] vector[K] mean_exposure; // means per category & group
  array[M, L] vector[K] css_exposure; // centered sum of squares per cat & group
  int n_test; // unique test locations & exposure groups
  array[n_test] vector[K] x_test; // stimulus locations in cue space
  array[n_test] int y_test; // group indicator of test trials
  array[n_test, M] int z_test_counts; // response counts per category M
  int<lower=0, upper=1> category_prior_known;
  array[category_prior_known ? 1 : 0] simplex[M] category_prior_data;
  int<lower=0, upper=1> weights_known; // 1 if group weights fixed
  array[weights_known ? L : 0] simplex[K] weights_data;
  int<lower=0, upper=1> lapse_rate_known;
  array[lapse_rate_known ? 1 : 0] real<lower=0, upper=1> lapse_rate_data;
  int<lower=0, upper=1> lapse_bias_known;
  array[lapse_bias_known ? 1 : 0] simplex[M] lapse_bias_data;
  // Category moments (partial fixing support)
  int<lower=0, upper=M> M_m_0_known;
  array[M_m_0_known] int<lower=1, upper=M> idx_m_0_known;
  array[M_m_0_known] vector[K] m_0_known_data;
  int<lower=0, upper=M> M_S_0_known;
  array[M_S_0_known] int<lower=1, upper=M> idx_S_0_known;
// 1 = S_0, 2 = E[Sigma_0]
  array[M_S_0_known] int<lower=1, upper=2> mode_S_0_known;
  array[M_S_0_known] cov_matrix[K] S_0_known_data;

  // Truly known category moments
  array[M] int<lower=0, upper=1> is_mu_0_known;
  array[M] vector[K] mu_0_known_data;

  array[M] int<lower=0, upper=1> is_Sigma_0_known;
  array[M] cov_matrix[K] Sigma_0_known_data;
  // Pseudocounts (shared vs. category-specific and partial fixing)
  int<lower=0, upper=1> kappa_0_category_specific;
  int<lower=0, upper=M> M_kappa_0_known;
  array[M_kappa_0_known] int<lower=1, upper=M> idx_kappa_0_known;
  array[M_kappa_0_known] real<lower=K> kappa_0_known_data;
  int<lower=0, upper=1> nu_0_category_specific;
  int<lower=0, upper=M> M_nu_0_known;
  array[M_nu_0_known] int<lower=1, upper=M> idx_nu_0_known;
  array[M_nu_0_known] real<lower=K+1> nu_0_known_data;
  // Prior hyperparameters (standardized prior_{hyperparameter}_of_{par})
  vector<lower=0>[M] prior_alphas_of_category_prior;
  vector<lower=0>[K] prior_alphas_of_weights;
  real<lower=0> prior_alpha_of_lapse_rate;
  real<lower=0> prior_beta_of_lapse_rate;
  vector<lower=0>[M] prior_alphas_of_lapse_bias;
  array[M] vector[K] prior_mu_of_m;
  vector[K] prior_loc_of_tau;
  vector<lower=0>[K] prior_scale_of_tau;
  array[kappa_0_category_specific ? M : 1] real prior_loc_of_kappa;
  array[kappa_0_category_specific ? M : 1] real<lower=0> prior_scale_of_kappa;
  array[nu_0_category_specific ? M : 1] real prior_loc_of_nu;
  array[nu_0_category_specific ? M : 1] real<lower=0> prior_scale_of_nu;
  matrix[K, K] INV_SCALE; // inverse affine transform matrix
  vector[K] shift; // affine transform center
  int<lower=0, upper=1> split_loglik_per_observation;
}

transformed data {
  // Compute unknown counts and indices for m_0
  int M_m_0_unknown = M - M_m_0_known;
  array[M_m_0_unknown] int<lower=1, upper=M> idx_m_0_unknown;
  {
    int count = 0;
    for (m in 1:M) {
      int is_known = 0;
      for (k in 1:M_m_0_known) {
        if (idx_m_0_known[k] == m) {
          is_known = 1;
          break;
        }
      }
      if (!is_known) {
        count += 1;
        idx_m_0_unknown[count] = m;
      }
    }
  }
  // Compute unknown counts and indices for S_0
  int M_S_0_unknown = M - M_S_0_known;
  array[M_S_0_unknown] int<lower=1, upper=M> idx_S_0_unknown;
  {
    int count = 0;
    for (m in 1:M) {
      int is_known = 0;
      for (k in 1:M_S_0_known) {
        if (idx_S_0_known[k] == m) {
          is_known = 1;
          break;
        }
      }
      if (!is_known) {
        count += 1;
        idx_S_0_unknown[count] = m;
      }
    }
  }
  // Compute unknown counts and indices for kappa_0
  int M_kappa_0_unknown;
  if (kappa_0_category_specific) {
    M_kappa_0_unknown = M - M_kappa_0_known;
  } else {
    M_kappa_0_unknown = (M_kappa_0_known > 0) ? 0 : 1;
  }
  array[kappa_0_category_specific ? M_kappa_0_unknown : 0]
    int<lower=1, upper=M> idx_kappa_0_unknown;
  if (kappa_0_category_specific) {
    int count = 0;
    for (m in 1:M) {
      int is_known = 0;
      for (k in 1:M_kappa_0_known) {
        if (idx_kappa_0_known[k] == m) {
          is_known = 1;
          break;
        }
      }
      if (!is_known) {
        count += 1;
        idx_kappa_0_unknown[count] = m;
      }
    }
  }
  // Compute unknown counts and indices for nu_0
  int M_nu_0_unknown;
  if (nu_0_category_specific) {
    M_nu_0_unknown = M - M_nu_0_known;
  } else {
    M_nu_0_unknown = (M_nu_0_known > 0) ? 0 : 1;
  }
  array[nu_0_category_specific ? M_nu_0_unknown : 0]
    int<lower=1, upper=M> idx_nu_0_unknown;
  if (nu_0_category_specific) {
    int count = 0;
    for (m in 1:M) {
      int is_known = 0;
      for (k in 1:M_nu_0_known) {
        if (idx_nu_0_known[k] == m) {
          is_known = 1;
          break;
        }
      }
      if (!is_known) {
        count += 1;
        idx_nu_0_unknown[count] = m;
      }
    }
  }
}

parameters {
  array[weights_known ? 0 : L] simplex[K] weights_param;
  // Pseudocounts
  array[kappa_0_category_specific ? 0 : M_kappa_0_unknown]
    real<lower=K> kappa_0_shared_param;
  array[kappa_0_category_specific ? M_kappa_0_unknown : 0]
    real<lower=K> kappa_0_cat_param;
  array[nu_0_category_specific ? 0 : M_nu_0_unknown]
    real<lower=K+1> nu_0_shared_param;
  array[nu_0_category_specific ? M_nu_0_unknown : 0]
    real<lower=K+1> nu_0_cat_param;
  // Category probabilities and lapse
  array[category_prior_known ? 0 : 1] simplex[M] category_prior_param;
  array[lapse_rate_known ? 0 : 1] real<lower=0, upper=1> lapse_rate_param;
  array[lapse_bias_known ? 0 : 1] simplex[M] lapse_bias_param;
  // Category moments (sampled only for unknown categories)
  array[M_m_0_unknown] vector[K] m_0_unknown_param;
  vector<lower=0>[M_m_0_unknown > 0 ? K : 0] m_0_tau_transformed;
  array[M_S_0_unknown] vector<lower=0>[K] tau_0_param_transformed;
}

transformed parameters {
  array[L] simplex[K] weights = weights_known ? weights_data : weights_param;
  simplex[M] category_prior = category_prior_known ?
    category_prior_data[1] : category_prior_param[1];
  real lapse_rate = lapse_rate_known ?
    lapse_rate_data[1] : lapse_rate_param[1];
  simplex[M] lapse_bias = lapse_bias_known ?
    lapse_bias_data[1] : lapse_bias_param[1];
  // Reassemble kappa_0 and nu_0 for all M categories
  array[M] real<lower=K> kappa_0;
  if (kappa_0_category_specific) {
    for (i in 1:M_kappa_0_known) {
      kappa_0[idx_kappa_0_known[i]] = kappa_0_known_data[i];
    }
    for (i in 1:M_kappa_0_unknown) {
      kappa_0[idx_kappa_0_unknown[i]] = kappa_0_cat_param[i];
    }
  } else {
    real val = (M_kappa_0_known > 0) ?
      kappa_0_known_data[1] : kappa_0_shared_param[1];
    for (m in 1:M) kappa_0[m] = val;
  }
  array[M] real<lower=K+1> nu_0;
  if (nu_0_category_specific) {
    for (i in 1:M_nu_0_known) {
      nu_0[idx_nu_0_known[i]] = nu_0_known_data[i];
    }
    for (i in 1:M_nu_0_unknown) {
      nu_0[idx_nu_0_unknown[i]] = nu_0_cat_param[i];
    }
  } else {
    real val = (M_nu_0_known > 0) ?
      nu_0_known_data[1] : nu_0_shared_param[1];
    for (m in 1:M) nu_0[m] = val;
  }
  // Reassemble m_0_transformed
  array[M] vector[K] m_0_transformed;
  for (i in 1:M_m_0_known) {
    m_0_transformed[idx_m_0_known[i]] = m_0_known_data[i];
  }
  for (i in 1:M_m_0_unknown) {
    m_0_transformed[idx_m_0_unknown[i]] = m_0_unknown_param[i];
  }
  // Reassemble tau_0_transformed
  array[M] vector<lower=0>[K] tau_0_transformed;
  for (i in 1:M_S_0_known) {
    int cat = idx_S_0_known[i];
    if (mode_S_0_known[i] == 1) {
      tau_0_transformed[cat] = sqrt(diagonal(S_0_known_data[i]) / nu_0[cat]);
    } else {
      tau_0_transformed[cat] = sqrt(
        diagonal(S_0_known_data[i]) * ((nu_0[cat] - 2) / nu_0[cat])
      );
    }
  }
  for (i in 1:M_S_0_unknown) {
    int cat = idx_S_0_unknown[i];
    tau_0_transformed[cat] = tau_0_param_transformed[i];
  }
  // updated beliefs depend on input and group
  array[M, L] real<lower=K> kappa_n;
  array[M, L] real<lower=K+1> nu_n;
  array[M, L] vector[K] m_n_transformed;
  array[M, L] vector[K] tau_n_transformed;
  array[M, L] vector[K] t_scale_transformed;
  array[n_test] simplex[M] p_test_conj;
  array[n_test] vector[M] log_p_test_conj;

  for (cat in 1:M) {
    for (group in 1:L) {
      if (is_mu_0_known[cat] == 1 && is_Sigma_0_known[cat] == 1) {
        // Case 4: Known mu, Known Sigma
        m_n_transformed[cat, group] = mu_0_known_data[cat];
        tau_n_transformed[cat, group] = sqrt(
          diagonal(Sigma_0_known_data[cat])
        );
        t_scale_transformed[cat, group] = tau_n_transformed[cat, group];
      } else if (is_mu_0_known[cat] == 1) {
        // Case 3: Known mu, Unknown Sigma
        nu_n[cat, group] = nu_0[cat] + n_exposure[cat, group];
        m_n_transformed[cat, group] = mu_0_known_data[cat];
        tau_n_transformed[cat, group] = nix_known_mu_update_tau(
          nu_0[cat],
          tau_0_transformed[cat],
          n_exposure[cat, group],
          mean_exposure[cat, group],
          css_exposure[cat, group],
          mu_0_known_data[cat]
        );
        t_scale_transformed[cat, group] = nix_known_mu_predictive_scale(
          tau_n_transformed[cat, group]
        );
      } else if (is_Sigma_0_known[cat] == 1) {
        // Case 2: Unknown mu, Known Sigma
        kappa_n[cat, group] = kappa_0[cat] + n_exposure[cat, group];
        m_n_transformed[cat, group] = nix_update_mean(
          kappa_0[cat],
          m_0_transformed[cat],
          n_exposure[cat, group],
          mean_exposure[cat, group]
        );
        tau_n_transformed[cat, group] = sqrt(
          diagonal(Sigma_0_known_data[cat])
        );
        t_scale_transformed[cat, group] = nix_known_sigma_predictive_scale(
          tau_n_transformed[cat, group],
          kappa_n[cat, group]
        );
      } else {
        // Case 1: Unknown mu, Unknown Sigma (Standard NIX)
        kappa_n[cat, group] = kappa_0[cat] + n_exposure[cat, group];
        nu_n[cat, group] = nu_0[cat] + n_exposure[cat, group];
        m_n_transformed[cat, group] = nix_update_mean(
          kappa_0[cat],
          m_0_transformed[cat],
          n_exposure[cat, group],
          mean_exposure[cat, group]
        );
        tau_n_transformed[cat, group] = nix_update_tau(
          kappa_0[cat],
          nu_0[cat],
          tau_0_transformed[cat],
          n_exposure[cat, group],
          mean_exposure[cat, group],
          css_exposure[cat, group],
          m_0_transformed[cat]
        );
        t_scale_transformed[cat, group] = nix_predictive_scale(
          tau_n_transformed[cat, group],
          kappa_n[cat, group]
        );
      }
    }
  }
  
  // compute posterior probabilities for each test stimulus (group-weighted)
  for (j in 1:n_test) {
    int group = y_test[j];
    for (cat in 1:M) {
      vector[K] cue_log_lik;
      if (is_Sigma_0_known[cat] == 1) {
        for (k in 1:K) {
          cue_log_lik[k] = normal_lpdf(
            x_test[j][k] |
            m_n_transformed[cat, group][k],
            t_scale_transformed[cat, group][k]
          );
        }
      } else {
        for (k in 1:K) {
          cue_log_lik[k] = student_t_lpdf(
            x_test[j][k] |
            nu_n[cat, group],
            m_n_transformed[cat, group][k],
            t_scale_transformed[cat, group][k]
          );
        }
      }
      log_p_test_conj[j, cat] = dot_product(weights[group], cue_log_lik) +
        log(category_prior[cat]);
    }
    p_test_conj[j] = exp(log_p_test_conj[j] - log_sum_exp(log_p_test_conj[j]));
  }
}

model {
  if (!weights_known) {
    for (g in 1:L) {
      weights_param[g] ~ dirichlet(prior_alphas_of_weights);
    }
  }
  // Pseudocount priors
  if (kappa_0_category_specific) {
    for (i in 1:M_kappa_0_unknown) {
      int cat = idx_kappa_0_unknown[i];
      kappa_0_cat_param[i] ~ cauchy(
        prior_loc_of_kappa[cat], prior_scale_of_kappa[cat]
      );
    }
  } else {
    if (M_kappa_0_unknown > 0) {
      kappa_0_shared_param[1] ~ cauchy(
        prior_loc_of_kappa[1], prior_scale_of_kappa[1]
      );
    }
  }
  if (nu_0_category_specific) {
    for (i in 1:M_nu_0_unknown) {
      int cat = idx_nu_0_unknown[i];
      nu_0_cat_param[i] ~ cauchy(
        prior_loc_of_nu[cat], prior_scale_of_nu[cat]
      );
    }
  } else {
    if (M_nu_0_unknown > 0) {
      nu_0_shared_param[1] ~ cauchy(
        prior_loc_of_nu[1], prior_scale_of_nu[1]
      );
    }
  }
  // Category probabilities and lapse
  if (!category_prior_known) {
    category_prior_param[1] ~ dirichlet(prior_alphas_of_category_prior);
  }
  if (!lapse_rate_known) {
    lapse_rate_param[1] ~ beta(
      prior_alpha_of_lapse_rate, prior_beta_of_lapse_rate
    );
  }
  if (!lapse_bias_known) {
    lapse_bias_param[1] ~ dirichlet(prior_alphas_of_lapse_bias);
  }
  // Prior for m_0
  if (M_m_0_unknown > 0) {
    m_0_tau_transformed ~ cauchy(prior_loc_of_tau, prior_scale_of_tau);
    for (i in 1:M_m_0_unknown) {
      int cat = idx_m_0_unknown[i];
      m_0_unknown_param[i] ~ normal(
        prior_mu_of_m[cat], m_0_tau_transformed
      );
    }
  }
  // Prior for tau_0
  if (M_S_0_unknown > 0) {
    for (i in 1:M_S_0_unknown) {
      tau_0_param_transformed[i] ~ cauchy(
        prior_loc_of_tau, prior_scale_of_tau
      );
    }
  }
  // Likelihood
  for (i in 1:n_test) {
    z_test_counts[i] ~ multinomial(
      lapse_mixture(p_test_conj[i], lapse_rate, lapse_bias)
    );
  }
}

generated quantities {
#include include/gq_multinomial_loglik.stan
  // Backtransformed parameters in natural physical units (clean base names)
  array[M] vector[K] m_0;
  array[M] vector[K] S_0;
  array[M, L] vector[K] m_n;
  array[M, L] vector[K] S_n;
  for (cat in 1:M) {
    m_0[cat] = unscale_mean(INV_SCALE, shift, m_0_transformed[cat]);
    S_0[cat] = unscale_diag_var(INV_SCALE, tau_0_transformed[cat]);
    for (group in 1:L) {
      m_n[cat, group] = unscale_mean(
        INV_SCALE, shift, m_n_transformed[cat, group]
      );
      S_n[cat, group] = unscale_diag_var(
        INV_SCALE, tau_n_transformed[cat, group]
      );
    }
  }
}
