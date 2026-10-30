/*
  Weighted Integrator - MNIX Ideal Adaptor

  Fit multinomial response data using a belief-updating model to infer prior
  parameters and global projection weights under early cue integration.

  Cognitive & Biological Interpretation:
  Under early integration with a weighted integrator, cues are projected into
  a single composite dimension:
    z = w' x = sum_{k=1}^K w_k * x_k
  followed by 1D Normal-Inverse-Chisquare (NIX) adaptation and
  categorization on z. The global projection weights w in Delta^{K-1} reflect an
  attention-weighted or sensorily filtered linear combination of cues shared
  across categories and groups.

  Backtransformation of Parameters:
  In generated quantities, the projection weights are transformed back to the
  physical cue space gradient using the affine transformation matrix:
    w_orig = SCALE' * weights_transformed
    weights = w_orig / sum(w_orig)
  and the composite category parameters on z are backtransformed via:
    s_w = sqrt(w' * INV_SCALE * INV_SCALE' * w)
    shift_w = dot_product(w, shift)
    m_0 = s_w * m_0_transformed - shift_w
    S_0 = s_w^2 * S_0_transformed
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
  array[M, L] matrix[K, K] css_exposure; // scatter matrix of centered squares
  int n_test; // unique test locations & exposure groups
  array[n_test] vector[K] x_test; // stimulus locations in cue space
  array[n_test] int y_test; // group indicator of test trials
  array[n_test, M] int z_test_counts; // response counts per category M
  int<lower=0, upper=1> category_prior_known;
  array[category_prior_known ? 1 : 0] simplex[M] category_prior_data;
  int<lower=0, upper=1> weights_known; // 1 if global weights fixed
  array[weights_known ? 1 : 0] simplex[K] weights_data;
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
  array[M_kappa_0_known] real<lower=1> kappa_0_known_data;
  int<lower=0, upper=1> nu_0_category_specific;
  int<lower=0, upper=M> M_nu_0_known;
  array[M_nu_0_known] int<lower=1, upper=M> idx_nu_0_known;
  array[M_nu_0_known] real<lower=2> nu_0_known_data;
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
  matrix[K, K] SCALE = inverse(INV_SCALE);
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
  array[weights_known ? 0 : 1] simplex[K] weights_param_transformed;
  // Pseudocounts
  array[kappa_0_category_specific ? 0 : M_kappa_0_unknown]
    real<lower=1> kappa_0_shared_param;
  array[kappa_0_category_specific ? M_kappa_0_unknown : 0]
    real<lower=1> kappa_0_cat_param;
  array[nu_0_category_specific ? 0 : M_nu_0_unknown]
    real<lower=2> nu_0_shared_param;
  array[nu_0_category_specific ? M_nu_0_unknown : 0]
    real<lower=2> nu_0_cat_param;
  // Category probabilities and lapse
  array[category_prior_known ? 0 : 1] simplex[M] category_prior_param;
  array[lapse_rate_known ? 0 : 1] real<lower=0, upper=1> lapse_rate_param;
  array[lapse_bias_known ? 0 : 1] simplex[M] lapse_bias_param;
  // Category moments on composite z (sampled only for unknown categories)
  array[M_m_0_unknown] real m_0_unknown_param;
  vector<lower=0>[M_m_0_unknown > 0 ? 1 : 0] m_0_tau_transformed;
  array[M_S_0_unknown] real<lower=0> tau_0_param_transformed;
}

transformed parameters {
  simplex[K] weights_transformed = weights_known ?
    weights_data[1] : weights_param_transformed[1];
  simplex[M] category_prior = category_prior_known ?
    category_prior_data[1] : category_prior_param[1];
  real lapse_rate = lapse_rate_known ?
    lapse_rate_data[1] : lapse_rate_param[1];
  simplex[M] lapse_bias = lapse_bias_known ?
    lapse_bias_data[1] : lapse_bias_param[1];
  real prior_loc_of_tau_1d = dot_product(
    weights_transformed, prior_loc_of_tau
  );
  real prior_scale_of_tau_1d = dot_product(
    weights_transformed, prior_scale_of_tau
  );
  array[M] real prior_mu_of_m_1d;
  for (cat in 1:M) {
    prior_mu_of_m_1d[cat] = dot_product(
      weights_transformed, prior_mu_of_m[cat]
    );
  }
  // Precomputed 1D composite statistics on z = w' x
  array[M, L] real z_mean_exposure;
  array[M, L] real z_css_exposure;
  for (cat in 1:M) {
    for (group in 1:L) {
      z_mean_exposure[cat, group] = dot_product(
        weights_transformed, mean_exposure[cat, group]
      );
      z_css_exposure[cat, group] = quad_form(
        css_exposure[cat, group], weights_transformed
      );
    }
  }
  array[n_test] real z_test;
  for (j in 1:n_test) {
    z_test[j] = dot_product(weights_transformed, x_test[j]);
  }
  // Reassemble kappa_0 and nu_0 for all M categories
  array[M] real<lower=1> kappa_0;
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
  array[M] real<lower=2> nu_0;
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
  // Reassemble m_0_transformed on composite z
  array[M] real m_0_transformed;
  for (i in 1:M_m_0_known) {
    m_0_transformed[idx_m_0_known[i]] = dot_product(
      weights_transformed, m_0_known_data[i]
    );
  }
  for (i in 1:M_m_0_unknown) {
    m_0_transformed[idx_m_0_unknown[i]] = m_0_unknown_param[i];
  }
  // Reassemble S_0_transformed on composite z
  array[M] real<lower=0> S_0_transformed;
  for (i in 1:M_S_0_known) {
    int cat = idx_S_0_known[i];
    real S_1d = quad_form(S_0_known_data[i], weights_transformed);
    if (mode_S_0_known[i] == 1) {
      S_0_transformed[cat] = S_1d;
    } else {
      S_0_transformed[cat] = S_1d * (nu_0[cat] - 2);
    }
  }
  for (i in 1:M_S_0_unknown) {
    int cat = idx_S_0_unknown[i];
    S_0_transformed[cat] = tau_0_param_transformed[i]^2;
  }
  // updated 1D beliefs on z depend on input and group
  array[M, L] real<lower=1> kappa_n;
  array[M, L] real<lower=2> nu_n;
  array[M, L] real m_n_transformed;
  array[M, L] real<lower=0> S_n_transformed;
  array[M, L] real<lower=0> t_scale_transformed;
  array[n_test] simplex[M] p_test_conj;
  array[n_test] vector[M] log_p_test_conj;

  for (cat in 1:M) {
    for (group in 1:L) {
      if (is_mu_0_known[cat] == 1 && is_Sigma_0_known[cat] == 1) {
        // Case 4: Known mu, Known Sigma
        m_n_transformed[cat, group] = m_0_transformed[cat];
        S_n_transformed[cat, group] = S_0_transformed[cat];
        t_scale_transformed[cat, group] = sqrt(S_0_transformed[cat]);
      } else if (is_mu_0_known[cat] == 1) {
        // Case 3: Known mu, Unknown Sigma
        nu_n[cat, group] = nu_0[cat] + n_exposure[cat, group];
        m_n_transformed[cat, group] = m_0_transformed[cat];
        S_n_transformed[cat, group] = nix_1d_known_mu_update_S(
          S_0_transformed[cat],
          n_exposure[cat, group],
          z_mean_exposure[cat, group],
          z_css_exposure[cat, group],
          m_0_transformed[cat]
        );
        t_scale_transformed[cat, group] = nix_1d_known_mu_predictive_scale(
          S_n_transformed[cat, group],
          nu_n[cat, group]
        );
      } else if (is_Sigma_0_known[cat] == 1) {
        // Case 2: Unknown mu, Known Sigma
        kappa_n[cat, group] = kappa_0[cat] + n_exposure[cat, group];
        m_n_transformed[cat, group] = nix_1d_update_mean(
          kappa_0[cat],
          m_0_transformed[cat],
          n_exposure[cat, group],
          z_mean_exposure[cat, group]
        );
        S_n_transformed[cat, group] = S_0_transformed[cat];
        t_scale_transformed[cat, group] = nix_1d_known_sigma_predictive_scale(
          sqrt(S_0_transformed[cat]),
          kappa_n[cat, group]
        );
      } else {
        // Case 1: Unknown mu, Unknown Sigma (Standard 1D NIX)
        kappa_n[cat, group] = kappa_0[cat] + n_exposure[cat, group];
        nu_n[cat, group] = nu_0[cat] + n_exposure[cat, group];
        m_n_transformed[cat, group] = nix_1d_update_mean(
          kappa_0[cat],
          m_0_transformed[cat],
          n_exposure[cat, group],
          z_mean_exposure[cat, group]
        );
        S_n_transformed[cat, group] = nix_1d_update_S(
          kappa_0[cat],
          S_0_transformed[cat],
          n_exposure[cat, group],
          z_mean_exposure[cat, group],
          z_css_exposure[cat, group],
          m_0_transformed[cat]
        );
        t_scale_transformed[cat, group] = nix_1d_predictive_scale(
          S_n_transformed[cat, group],
          kappa_n[cat, group],
          nu_n[cat, group]
        );
      }
    }
  }
  
  // compute posterior probabilities for each test stimulus on composite z
  for (j in 1:n_test) {
    int group = y_test[j];
    for (cat in 1:M) {
      if (is_Sigma_0_known[cat] == 1) {
        log_p_test_conj[j, cat] = normal_lpdf(
          z_test[j] |
          m_n_transformed[cat, group],
          t_scale_transformed[cat, group]
        ) + log(category_prior[cat]);
      } else {
        log_p_test_conj[j, cat] = student_t_lpdf(
          z_test[j] |
          nu_n[cat, group],
          m_n_transformed[cat, group],
          t_scale_transformed[cat, group]
        ) + log(category_prior[cat]);
      }
    }
    p_test_conj[j] = exp(log_p_test_conj[j] - log_sum_exp(log_p_test_conj[j]));
  }
}

model {
  if (!weights_known) {
    weights_param_transformed[1] ~ dirichlet(prior_alphas_of_weights);
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
  // Prior for m_0 on composite z
  if (M_m_0_unknown > 0) {
    m_0_tau_transformed[1] ~ cauchy(
      prior_loc_of_tau_1d, prior_scale_of_tau_1d
    );
    for (i in 1:M_m_0_unknown) {
      int cat = idx_m_0_unknown[i];
      m_0_unknown_param[i] ~ normal(
        prior_mu_of_m_1d[cat], m_0_tau_transformed[1]
      );
    }
  }
  // Prior for tau_0 on composite z
  if (M_S_0_unknown > 0) {
    for (i in 1:M_S_0_unknown) {
      tau_0_param_transformed[i] ~ cauchy(
        prior_loc_of_tau_1d, prior_scale_of_tau_1d
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
  // Backtransformed parameters in natural physical units
  vector[K] w_orig = SCALE' * weights_transformed;
  simplex[K] weights = w_orig / sum(w_orig);
  real s_w = sqrt(quad_form(INV_SCALE * INV_SCALE', weights));
  real shift_w = dot_product(weights, shift);
  array[M] real m_0;
  array[M] real S_0;
  array[M, L] real m_n;
  array[M, L] real S_n;
  for (cat in 1:M) {
    m_0[cat] = unscale_1d_mean(s_w, shift_w, m_0_transformed[cat]);
    S_0[cat] = unscale_1d_scatter(s_w, S_0_transformed[cat]);
    for (group in 1:L) {
      m_n[cat, group] = unscale_1d_mean(
        s_w, shift_w, m_n_transformed[cat, group]
      );
      S_n[cat, group] = unscale_1d_scatter(
        s_w, S_n_transformed[cat, group]
      );
    }
  }
}
