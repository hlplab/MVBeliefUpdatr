/*
  Category-Weighted Integrator - MNIX Ideal Adaptor

  Fit multinomial response data using a belief-updating model to infer prior
  parameters and category-specific projection weights under early cue integration.

  Cognitive & Biological Interpretation:
  Under early integration with a category-weighted integrator, each category c projects
  multi-cue stimuli into its own composite dimension:
    z_c = w_c' x = sum_{k=1}^K w_{c, k} * x_k
  followed by 1D Normal-Inverse-Chisquare (NIX) adaptation on z_c.
  This reflects perceptual mechanisms where category templates filter inputs along
  category-specific diagnostic axes (e.g. phoneme-specific phonetic projections).

  Backtransformation of Parameters:
  In generated quantities, category projection weights are transformed back to
  physical cue gradients:
    w_{c, orig} = SCALE' * weights_transformed[c]
    weights[c] = w_{c, orig} / sum(w_{c, orig})
  and composite parameters on z_c are backtransformed via category-specific scalings:
    s_{w, c} = sqrt(w_c' * INV_SCALE * INV_SCALE' * w_c)
    shift_{w, c} = dot_product(w_c, shift)
    m_0[c] = s_{w, c} * m_0_transformed[c] - shift_{w, c}
    S_0[c] = s_{w, c}^2 * S_0_transformed[c]
*/

functions {
#include include/mvbu_functions.stan
}


data {
  int K;                                                   // number of features (cues >= 2)
  int L;                                                   // number of exposure groups (e.g. subjects)
  int M;                                                   // number of categories

  array[M,L] int<lower=0> N_exposure;                      // observations per category (M) and group (L)
  array[M,L] vector[K] x_mean_exposure;                    // means per category (M) and group (L)
  array[M,L] matrix[K,K] x_css_exposure;                   // scatter matrix of centered squares

  int N_test;                                              // unique test locations & exposure groups
  array[N_test] vector[K] x_test;                          // stimulus locations in cue space
  array[N_test] int y_test;                                // group indicator of test trials
  array[N_test,M] int z_test_counts;                       // response counts per category M

  simplex[M] prior_p;                                      // prior category probabilities

  int<lower=0, upper=1> weights_known;                     // 1 if weights fixed by user, 0 if estimated
  int<lower=0, upper=1> lapse_rate_known;
  int<lower=0, upper=1> lapse_bias_known;
  int<lower=0, upper=1> mu_0_known;
  int<lower=0, upper=1> Sigma_0_known;
  array[weights_known ? M : 0] simplex[K] weights_data;    // optional: user-provided category weights
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
  matrix[K, K] SCALE = inverse(INV_SCALE);
  simplex[M] lapse_bias = lapse_bias_known ? lapse_bias_data[1] : rep_vector(1.0 / M, M);
  real<lower=0> sigma_kappanu = max(max(to_array_1d(N_exposure)) * 4, 10);
  real tau_scale_1d = mean(tau_scale);
}

parameters {
  array[weights_known ? 0 : M] simplex[K] weights_param_transformed; // category projection weights (standardized)
  real<lower=1> kappa_0;                                             // prior pseudocount for composite mu
  real<lower=2> nu_0;                                                // prior pseudocount for composite Sigma

  array[lapse_rate_known ? 0 : 1] real<lower=0, upper=1> lapse_rate_param;

  array[mu_0_known ? 0 : M] real m_0_param_transformed;              // prior expected mean on z (standardized)
  array[mu_0_known ? 0 : 1] real<lower=0> m_0_tau_transformed;       // prior SD of m_0
  array[Sigma_0_known ? 0 : M] real<lower=0> tau_0_param_transformed; // prior scale on z (standardized)
}

transformed parameters {
  array[M] simplex[K] weights_transformed = weights_known ? weights_data : weights_param_transformed;

  // Early-integrated 1D exposure statistics and test stimuli computed with category weights
  array[M,L] real z_mean_exposure;
  array[M,L] real z_css_exposure;
  for (cat in 1:M) {
    for (group in 1:L) {
      z_mean_exposure[cat,group] = dot_product(weights_transformed[cat], x_mean_exposure[cat,group]);
      z_css_exposure[cat,group] = quad_form(x_css_exposure[cat,group], weights_transformed[cat]);
    }
  }

  array[N_test, M] real z_test;
  for (j in 1:N_test) {
    for (cat in 1:M) {
      z_test[j, cat] = dot_product(weights_transformed[cat], x_test[j]);
    }
  }

  array[mu_0_known ? M : 0] real mu_0_1d_data;
  if (mu_0_known) {
    for (cat in 1:M) {
      mu_0_1d_data[cat] = dot_product(weights_transformed[cat], mu_0_data[cat]);
    }
  }

  array[Sigma_0_known ? M : 0] real Sigma_0_1d_data;
  if (Sigma_0_known) {
    for (cat in 1:M) {
      Sigma_0_1d_data[cat] = quad_form(Sigma_0_data[cat], weights_transformed[cat]);
    }
  }

  array[M] real m_0_transformed = mu_0_known ? mu_0_1d_data : m_0_param_transformed;
  array[M] real S_0_transformed;
  real lapse_rate = lapse_rate_known ? lapse_rate_data[1] : lapse_rate_param[1];

  // updated 1D beliefs on z depend on input and group
  array[M,L] real<lower=1> kappa_n;
  array[M,L] real<lower=2> nu_n;
  array[M,L] real m_n_transformed;
  array[M,L] real<lower=0> S_n_transformed;
  array[M,L] real<lower=0> t_scale_transformed;

  array[N_test] simplex[M] p_test_conj;
  array[N_test] vector[M] log_p_test_conj;

  for (cat in 1:M) {
    S_0_transformed[cat] = Sigma_0_known ? Sigma_0_1d_data[cat] * (nu_0 - 2) : tau_0_param_transformed[cat]^2;

    for (group in 1:L) {
      kappa_n[cat,group] = kappa_0 + N_exposure[cat,group];
      nu_n[cat,group] = nu_0 + N_exposure[cat,group];
      m_n_transformed[cat,group] =
        nix_1d_update_mean(kappa_0, m_0_transformed[cat], N_exposure[cat,group], z_mean_exposure[cat,group]);
      S_n_transformed[cat,group] =
        nix_1d_update_S(kappa_0, S_0_transformed[cat], N_exposure[cat,group], z_mean_exposure[cat,group], z_css_exposure[cat,group], m_0_transformed[cat]);
      t_scale_transformed[cat,group] =
        nix_1d_predictive_scale(S_n_transformed[cat,group], kappa_n[cat,group], nu_n[cat,group]);
    }
  }

  // compute posterior probabilities for each test stimulus on category composite dimensions
  for (j in 1:N_test) {
    int group = y_test[j];
    for (cat in 1:M) {
      log_p_test_conj[j,cat] =
        student_t_lpdf(z_test[j, cat] |
                       nu_n[cat,group],
                       m_n_transformed[cat,group],
                       t_scale_transformed[cat,group]) + log(prior_p[cat]);
    }
    p_test_conj[j] = exp(log_p_test_conj[j] - log_sum_exp(log_p_test_conj[j]));
  }
}

model {
  if (!weights_known) {
    for (cat in 1:M) {
      weights_param_transformed[cat] ~ dirichlet(rep_vector(1.0, K));
    }
  }
  kappa_0 ~ normal(0, sigma_kappanu);
  nu_0 ~ normal(0, sigma_kappanu);

  if (!mu_0_known) {
    m_0_tau_transformed[1] ~ cauchy(0, tau_scale_1d);
    m_0_param_transformed ~ normal(0, m_0_tau_transformed[1]);
  }

  if (!Sigma_0_known) {
    for (cat in 1:M) {
      tau_0_param_transformed[cat] ~ cauchy(0, tau_scale_1d);
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
  array[M] simplex[K] weights;

  for (cat in 1:M) {
    vector[K] w_orig = SCALE' * weights_transformed[cat];
    weights[cat] = w_orig / sum(w_orig);

    real s_w = sqrt(quad_form(INV_SCALE * INV_SCALE', weights_transformed[cat]));
    real shift_w = dot_product(weights_transformed[cat], shift);

    m_0[cat] = unscale_1d_mean(s_w, shift_w, m_0_transformed[cat]);
    S_0[cat] = unscale_1d_scatter(s_w, S_0_transformed[cat]);
    for (group in 1:L) {
      m_n[cat,group] = unscale_1d_mean(s_w, shift_w, m_n_transformed[cat,group]);
      S_n[cat,group] = unscale_1d_scatter(s_w, S_n_transformed[cat,group]);
    }
  }
}
