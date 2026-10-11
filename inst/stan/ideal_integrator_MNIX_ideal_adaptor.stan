/*
  Ideal Integrator - MNIX Ideal Adaptor

  Fit multinomial response data using a belief-updating model under normative
  early cue integration.

  Cognitive & Biological Interpretation:
  Under early integration with an ideal integrator, cues are combined equally and
  linearly into a composite cue:
    z = sum_{k=1}^K (1/K) * x_k
  followed by 1D Normal-Inverse-Chisquare (NIX) adaptation and Student-t categorization
  on composite dimension z.
  Because integration weights are strictly fixed to w_k = 1/K, all projections
  are precomputed once in transformed data, providing maximum computational efficiency.

  Backtransformation of Parameters:
  In generated quantities, the composite parameters on z are backtransformed to
  physical units using the linear projection of the affine coordinate transform:
    s_w = sqrt(w' * INV_SCALE * INV_SCALE' * w)
    shift_w = dot_product(w, shift)
    m_0 = s_w * m_0_transformed - shift_w
    S_0 = s_w^2 * S_0_transformed
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
  simplex[K] weights_const = rep_vector(1.0 / K, K);
  simplex[M] lapse_bias = lapse_bias_known ? lapse_bias_data[1] : rep_vector(1.0 / M, M);

  // Precomputed 1D composite statistics on z = w' x
  array[M,L] real z_mean_exposure;
  array[M,L] real z_css_exposure;
  for (cat in 1:M) {
    for (group in 1:L) {
      z_mean_exposure[cat,group] = dot_product(weights_const, x_mean_exposure[cat,group]);
      z_css_exposure[cat,group] = quad_form(x_css_exposure[cat,group], weights_const);
    }
  }

  array[N_test] real z_test;
  for (j in 1:N_test) {
    z_test[j] = dot_product(weights_const, x_test[j]);
  }

  array[mu_0_known ? M : 0] real mu_0_1d_data;
  if (mu_0_known) {
    for (cat in 1:M) {
      mu_0_1d_data[cat] = dot_product(weights_const, mu_0_data[cat]);
    }
  }

  array[Sigma_0_known ? M : 0] real Sigma_0_1d_data;
  if (Sigma_0_known) {
    for (cat in 1:M) {
      Sigma_0_1d_data[cat] = quad_form(Sigma_0_data[cat], weights_const);
    }
  }

  real<lower=0> sigma_kappanu = max(max(to_array_1d(N_exposure)) * 4, 10);
  real tau_scale_1d = dot_product(weights_const, tau_scale);

  real s_w = sqrt(quad_form(INV_SCALE * INV_SCALE', weights_const));
  real shift_w = dot_product(weights_const, shift);
}

parameters {
  real<lower=1> kappa_0;                                        // prior pseudocount for composite mu
  real<lower=2> nu_0;                                           // prior pseudocount for composite Sigma

  array[lapse_rate_known ? 0 : 1] real<lower=0, upper=1> lapse_rate_param;

  array[mu_0_known ? 0 : M] real m_0_param_transformed;        // prior expected mean on z (standardized)
  array[mu_0_known ? 0 : 1] real<lower=0> m_0_tau_transformed; // prior SD of m_0
  array[Sigma_0_known ? 0 : M] real<lower=0> tau_0_param_transformed; // prior scale on z (standardized)
}

transformed parameters {
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

  // compute posterior probabilities for each test stimulus on composite dimension z
  for (j in 1:N_test) {
    int group = y_test[j];
    for (cat in 1:M) {
      log_p_test_conj[j,cat] =
        student_t_lpdf(z_test[j] |
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

  // Backtransformed parameters on composite cue in natural physical units (clean base names)
  array[M] real m_0;
  array[M] real S_0;
  array[M,L] real m_n;
  array[M,L] real S_n;
  simplex[K] weights = weights_const;

  for (cat in 1:M) {
    m_0[cat] = unscale_1d_mean(s_w, shift_w, m_0_transformed[cat]);
    S_0[cat] = unscale_1d_scatter(s_w, S_0_transformed[cat]);
    for (group in 1:L) {
      m_n[cat,group] = unscale_1d_mean(s_w, shift_w, m_n_transformed[cat,group]);
      S_n[cat,group] = unscale_1d_scatter(s_w, S_n_transformed[cat,group]);
    }
  }
}
