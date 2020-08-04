/*
 * Fit multinomial response data using a belief-updating model to infer prior
 * parameters.  A normal-inverse-Wishart prior is used, and it's assumed
 * that the subject knows the true labels of all the input stimuli (e.g.,
 * doesn't model any uncertainty in classification.
 *
 * This version has a lapse rate parameter (probability of random guessing)
 *
 * Input is in the form of raw data points for observations.
 *
 *
 * Modified by
 * Florian Jaeger
 * August 2020
 */

data {
  int M;                        // number of categories
  int L;                        // number of grouping levels (e.g. subjects)
  int K;                        // number of features

  /* Efficiency considerations: declaring n as matrix rather than array of int is supposedly faster
     (https://mc-stan.org/docs/2_18/stan-users-guide/indexing-efficiency-section.html). It also means
     that standard arithmetic functions can be directly applied to that matrix (unlike to more-than-1-
     dimensional arrays). */
  matrix[M,L] N;                // number of observations per category (m) and subject (l)
  vector[K] x_mean[M,L];        // means for each category (m) and subject (l)
  cov_matrix[K] x_ss[M,L];      // sum of uncentered squares matrix for each category (m) and subject (l)

  int N_test;                   // number of test trials
  vector[K] x_test[N_test];     // locations of test trials
  int y_test[N_test];           // group label of test trials
  int z_test_counts[N_test,M];  // responses for test trials

  /* For now, this script assumes that the observations (cue vectors) are centered. The prior
     mean of mu_0 is set to 0. Same for the prior location parameter for the cauchy prior over
     the variance of mu_0 */
  real<lower=0> tau_scale;      // scale of cauchy prior for variances of mu_0 (set to zero to ignore)
  real<lower=0> L_omega_scale;  // scale of LKJ prior for correlation of variance of mu_0 (set to zero to ignore)
}

transformed data {
  real sigma_kappanu;

  sigma_kappanu = max(N) > 0 ? max(N) * 4 : 10;
}

parameters {
  // these are all shared across subjects (same prior beliefs):
  real<lower=K> kappa_0;                // prior pseudocount for mean
  real<lower=K> nu_0;                   // prior pseudocount for sd

  vector[K] mu_0[M];                    // prior expected mean
  vector<lower=0>[K] mu_0_tau;          // prior variances of mu_0
  cholesky_factor_corr[K] mu_0_L_omega; // prior correlations of variances of mu_0

  vector<lower=0>[K] tau_0[M];          // prior variances of sigma_0
  cholesky_factor_corr[K] L_omega_0[M]; // prior correlations of variances of sigma_0

  real<lower=0, upper=1> lapse_rate;
}

transformed parameters {
  // updated beliefs depend on input/subject
  real<lower=K> kappa_n[M,L];           // updated mean pseudocount
  real<lower=K> nu_n[M,L];              // updated sd pseudocount
  vector[K] mu_n[M,L];                  // updated expected mean
  cov_matrix[K] sigma_0[M];             // prior expected covariance matrix
  cov_matrix[K] sigma_n[M,L];           // updated expected scatter matrix
  cov_matrix[K] t_scale[M,L];           // scale matrix of predictive t distribution

  simplex[M] p_test_conj[N_test];
  vector[M] log_p_test_conj[N_test];

  // update NIW parameters according to conjugate updating rules are taken from
  // Murphy (2007, p. 136)
  for (cat in 1:M) {
    // Get sigma_0 from its components: correlation matrix and vector of standard deviations
    sigma_0[cat] = quad_form_diag(multiply_lower_tri_self_transpose(L_omega_0[cat]), tau_0[cat]);
    for (subj in 1:L) {
      if (N[cat,subj] > 0 ) {
        kappa_n[cat,subj] = kappa_0 + N[cat,subj];
        nu_n[cat,subj] = nu_0 + N[cat,subj];
        mu_n[cat,subj] = (kappa_0 * mu_0[cat] + N[cat,subj] * x_mean[cat,subj]) /
                          kappa_n[cat,subj];
        sigma_n[cat,subj] = sigma_0[cat] +
                            x_ss[cat,subj] +
                            kappa_0 * mu_0[cat] * mu_0[cat]' -
                            kappa_n[cat,subj] * mu_n[cat,subj] * mu_n[cat,subj]';
      } else {
        kappa_n[cat,subj] = kappa_0;
        nu_n[cat,subj] = nu_0;
        mu_n[cat,subj] = mu_0[cat];
        sigma_n[cat,subj] = sigma_0[cat];
      }

      t_scale[cat,subj] = sigma_n[cat,subj] * (kappa_n[cat,subj] + 1) /
                                              (kappa_n[cat,subj] * (nu_n[cat,subj] - K + 1));
    }
  }

  // compute category probabilities for each of the test stimuli
  for (j in 1:N_test) {
    int subj;
    subj = y_test[j];
    // calculate un-normalized log prob for each category
    for (cat in 1:M) {
      log_p_test_conj[j,cat] = multi_student_t_lpdf(x_test[j] |
                                              nu_n[cat,subj] - K + 1,
                                              mu_n[cat,subj],
                                              t_scale[cat,subj]);
    }
    // normalize and store actual probs in simplex
    p_test_conj[j] = exp(log_p_test_conj[j] - log_sum_exp(log_p_test_conj[j]));
  }
}

model {
  vector[M] lapsing_probs;

  lapsing_probs = rep_vector(lapse_rate / M, M);

  /* Need to calculate category probabilities for each test trial. In order to deal with
     input that does not contain observations (in which case n_each == 0), we set the
     minimum value for SD to 10. */
  kappa_0 ~ normal(0, sigma_kappanu);
  nu_0 ~ normal(0, sigma_kappanu);

  /* Specifying prior for mu_0:
     - If no scale for variances (tau) of mu_0 is user-specified use weakly regularizing
       scale (5) for variances of mean.
     - If no scale for LKJ prior over correlation matrix of mu_0 is user-specified use
       scale 1 to set uniform prior over correlation matrices. */
  mu_0_tau ~ cauchy(0, tau_scale > 0 ? tau_scale : 5);
  mu_0_L_omega ~ lkj_corr_cholesky(L_omega_scale > 0 ? L_omega_scale : 1);
  mu_0 ~ multi_normal_cholesky(rep_vector(0, K), diag_pre_multiply(mu_0_tau, mu_0_L_omega));

  for (cat in 1:M) {
    tau_0[cat] ~ cauchy(0, tau_scale > 0 ? tau_scale : 10);
    L_omega_0[cat] ~ lkj_corr_cholesky(L_omega_scale > 0 ? L_omega_scale : 1);
  }

  for (i in 1:N_test) {
    z_test_counts[i] ~ multinomial(p_test_conj[i] * (1-lapse_rate) + lapsing_probs);
  }

}

generated quantities {
  matrix[K,K] mu_0_cor;
  matrix[K,K] mu_0_cov;

  mu_0_cor = multiply_lower_tri_self_transpose(mu_0_L_omega);
  mu_0_cov = quad_form_diag(mu_0_cor, mu_0_tau);
}
