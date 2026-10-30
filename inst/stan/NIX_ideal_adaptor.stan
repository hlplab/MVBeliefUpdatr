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
  int L; // number of groups (conditions or subjects)
  array[M, L] int<lower=0> n_exposure; // observations per cat & group
  matrix[M, L] mean_exposure;        // mean of training observations
  matrix[M, L] css_exposure; // sum of centered squares of observations
  int n_test;                         // number of test trials
  array[n_test] real x_test;          // locations of test trials
  array[n_test] int y_test;           // group labels for test trials
  array[n_test, M] int z_test_counts;  // responses for test trials

  // -------------------------------------------------------------------------
  // Parameters: Fixed data and flags
  // -------------------------------------------------------------------------
  int<lower=0, upper=1> category_prior_known;
  array[category_prior_known ? 1 : 0] simplex[M] category_prior_data;
  int<lower=0, upper=1> lapse_rate_known;
  array[lapse_rate_known ? 1 : 0] real<lower=0, upper=1> lapse_rate_data;
  int<lower=0, upper=1> lapse_bias_known;
  array[lapse_bias_known ? 1 : 0] simplex[M] lapse_bias_data;
  // Category moments (partial fixing support)
  int<lower=0, upper=M> M_m_0_known;
  array[M_m_0_known] int<lower=1, upper=M> idx_m_0_known;
  array[M_m_0_known] real m_0_known_data;
  int<lower=0, upper=M> M_S_0_known;
  array[M_S_0_known] int<lower=1, upper=M> idx_S_0_known;
// 1 = S_0, 2 = E[Sigma_0]
  array[M_S_0_known] int<lower=1, upper=2> mode_S_0_known;
  array[M_S_0_known] real<lower=0> S_0_known_data;
  // Pseudocounts (shared vs. category-specific and partial fixing)
  int<lower=0, upper=1> kappa_0_category_specific;
  int<lower=0, upper=M> M_kappa_0_known;
  array[M_kappa_0_known] int<lower=1, upper=M> idx_kappa_0_known;
  array[M_kappa_0_known] real<lower=1> kappa_0_known_data;
  int<lower=0, upper=1> nu_0_category_specific;
  int<lower=0, upper=M> M_nu_0_known;
  array[M_nu_0_known] int<lower=1, upper=M> idx_nu_0_known;
  array[M_nu_0_known] real<lower=2> nu_0_known_data;
  // -------------------------------------------------------------------------
  // Prior hyperparameters (standardized prior_{hyperparameter}_of_{par})
  // -------------------------------------------------------------------------
  vector<lower=0>[M] prior_alphas_of_category_prior;
  real<lower=0> prior_alpha_of_lapse_rate;
  real<lower=0> prior_beta_of_lapse_rate;
  vector<lower=0>[M] prior_alphas_of_lapse_bias;
  array[M] real prior_mu_of_m;
  real prior_loc_of_tau;
  real<lower=0> prior_scale_of_tau;
  array[kappa_0_category_specific ? M : 1] real prior_loc_of_kappa;
  array[kappa_0_category_specific ? M : 1] real<lower=0> prior_scale_of_kappa;
  array[nu_0_category_specific ? M : 1] real prior_loc_of_nu;
  array[nu_0_category_specific ? M : 1] real<lower=0> prior_scale_of_nu;
  // -------------------------------------------------------------------------
  // Transformation metadata
  // -------------------------------------------------------------------------
  real INV_SCALE;
  real shift;
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
  // Category moments (sampled only for unknown categories)
  array[M_m_0_unknown] real m_0_unknown_param;
  array[M_m_0_unknown > 0 ? 1 : 0] real<lower=0> m_0_tau_transformed;
  array[M_S_0_unknown] real<lower=0> tau_0_param_transformed;
}

transformed parameters {
  simplex[M] category_prior = category_prior_known ?
    category_prior_data[1] : category_prior_param[1];
  real lapse_rate = lapse_rate_known ?
    lapse_rate_data[1] : lapse_rate_param[1];
  simplex[M] lapse_bias = lapse_bias_known ?
    lapse_bias_data[1] : lapse_bias_param[1];
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
  // Reassemble m_0_transformed
  array[M] real m_0_transformed;
  for (i in 1:M_m_0_known) {
    m_0_transformed[idx_m_0_known[i]] = m_0_known_data[i];
  }
  for (i in 1:M_m_0_unknown) {
    m_0_transformed[idx_m_0_unknown[i]] = m_0_unknown_param[i];
  }
  // Reassemble S_0_transformed
  array[M] real<lower=0> S_0_transformed;
  for (i in 1:M_S_0_known) {
    int cat = idx_S_0_known[i];
    if (mode_S_0_known[i] == 1) {
      S_0_transformed[cat] = S_0_known_data[i];
    } else {
      S_0_transformed[cat] = S_0_known_data[i] * (nu_0[cat] - 2);
    }
  }
  for (i in 1:M_S_0_unknown) {
    int cat = idx_S_0_unknown[i];
    S_0_transformed[cat] = tau_0_param_transformed[i]^2;
  }
  // updated beliefs depend on input/group
  array[M, L] real<lower=1> kappa_n;
  array[M, L] real<lower=2> nu_n;
  array[M, L] real m_n_transformed;
  array[M, L] real<lower=0> S_n_transformed;
  array[M, L] real<lower=0> t_scale_transformed;
  array[n_test] simplex[M] p_test_conj;
  array[n_test] vector[M] log_p_test_conj;

  // update NIX2 parameters according to conjugate updating rules. Note that we're using the 
  // scatter parameterization (S_0, S_n) rather than the variance/scale parameterization of 
  // Murphy (2012). This keeps things parallel to the NIW model, and also means that post-
  // processing functions like get_expected_sigma_from_S have the same meaning for both 
  // types of models (rather than double transforming the NIX scale parameter when the goal 
  // is to extract the expected category variance).
  for (cat in 1:M) {
    for (group in 1:L) {
      kappa_n[cat, group] = kappa_0[cat] + n_exposure[cat, group];
      nu_n[cat, group] = nu_0[cat] + n_exposure[cat, group];
      m_n_transformed[cat, group] = nix_1d_update_mean(
        kappa_0[cat],
        m_0_transformed[cat],
        n_exposure[cat, group],
        mean_exposure[cat, group]
      );
      S_n_transformed[cat, group] = nix_1d_update_S(
        kappa_0[cat],
        S_0_transformed[cat],
        n_exposure[cat, group],
        mean_exposure[cat, group],
        css_exposure[cat, group],
        m_0_transformed[cat]
      );
      t_scale_transformed[cat, group] = nix_1d_predictive_scale(
        S_n_transformed[cat, group],
        kappa_n[cat, group],
        nu_n[cat, group]
      );
    }
  }

  // compute category probabilities for each test stimulus
  for (j in 1:n_test) {
    int group = y_test[j];
    for (cat in 1:M) {
      log_p_test_conj[j, cat] = student_t_lpdf(
        x_test[j] |
        nu_n[cat, group],
        m_n_transformed[cat, group],
        t_scale_transformed[cat, group]
      ) + log(category_prior[cat]);
    }
    p_test_conj[j] = exp(log_p_test_conj[j] - log_sum_exp(log_p_test_conj[j]));
  }
}

model {
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
    m_0_tau_transformed[1] ~ cauchy(prior_loc_of_tau, prior_scale_of_tau);
    for (i in 1:M_m_0_unknown) {
      int cat = idx_m_0_unknown[i];
      m_0_unknown_param[i] ~ normal(prior_mu_of_m[cat], m_0_tau_transformed[1]);
    }
  }
  // Prior for components of S_0
  if (M_S_0_unknown > 0) {
    for (i in 1:M_S_0_unknown) {
      tau_0_param_transformed[i] ~ cauchy(
        prior_loc_of_tau, prior_scale_of_tau
      );
    }
  }
  // Likelihood
  for (n in 1:n_test) {
    z_test_counts[n] ~ multinomial(
      lapse_mixture(p_test_conj[n], lapse_rate, lapse_bias)
    );
  }
}

generated quantities {
#include include/gq_multinomial_loglik.stan

  // Backtransformed parameters in natural physical units (clean base names)
  array[M] real m_0;
  array[M] real S_0;
  array[M, L] real m_n;
  array[M, L] real S_n;
  for (cat in 1:M) {
    m_0[cat] = unscale_1d_mean(INV_SCALE, shift, m_0_transformed[cat]);
    S_0[cat] = unscale_1d_scatter(INV_SCALE, S_0_transformed[cat]);
    for (group in 1:L) {
      m_n[cat, group] = unscale_1d_mean(
        INV_SCALE, shift, m_n_transformed[cat, group]
      );
      S_n[cat, group] = unscale_1d_scatter(
        INV_SCALE, S_n_transformed[cat, group]
      );
    }
  }
}
