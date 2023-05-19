/*
 * Fit multinomial response data using a belief-updating model to infer prior
 * parameters.  A normal-inverse-Wishart prior is used, and it's assumed
 * that the subject knows the true labels of all the input stimuli (e.g.,
 * doesn't model any uncertainty in classification.
 *
 * This version has a lapse rate parameter (probability of random guessing)
 *
 * Input is in the form of sufficient statistics of exposure observations, and
 * counts of responses for test observations.
 *
 *
 * Modified by Florian Jaeger, May 2023 to:
 *
 *    + Allow independent updating over each dimension of multivariate categories.
 */

data {
  int M;                        // number of categories
  int L;                        // number of grouping levels (e.g. subjects)
  int K;                        // number of features

  /* Efficiency considerations: declaring n as matrix rather than array of int is supposedly faster
     (https://mc-stan.org/docs/2_18/stan-users-guide/indexing-efficiency-section.html). It also means
     that standard arithmetic functions can be directly applied to that matrix (unlike to more-than-1-
     dimensional arrays). */
  matrix[M,L] N;                // number of observations per category (m) and group (l) in exposure
  vector[K] x_mean[M,L];        // means for each category (m) and group (l) in exposure
  vector[K] x_sd[M,L];          // sample standard deviation for each category (m) and group (l) in exposure

  int N_test;                   // number of test trials
  vector[K] x_test[N_test];     // locations of test trials
  int y_test[N_test];           // group label of test trials
  int z_test_counts[N_test,M];  // responses for test trials

  int<lower=0, upper=1> multiple_kappa_0;
  int<lower=0, upper=1> multiple_nu_0;

  int<lower=0, upper=1> m_0_known;
  int<lower=0, upper=1> S_0_known;
  vector[m_0_known ? K : 0] m_0_data[m_0_known ? M : 0];    // optional: user provided m_0 (prior mean of means)
  vector[S_0_known ? K : 0] S_0_data[S_0_known ? M : 0];    // optional: user provided S_0 (prior expected standard deviation)

  /* For now, this script assumes that the observations (cue vectors) are centered. The prior
     mean of m_0 is set to 0. Same for the prior location parameter for the cauchy prior over
     the variance of m_0 */
  // separate taus for each feature to capture that features can be on separate scales:
  // vector<lower=0>[K] tau_scales;  // scales of cauchy prior for variances along the K features (set to zero to ignore)
  real<lower=0> tau_scale;      // scale of cauchy prior for variances of m_0 (set to zero to ignore)
  real<lower=0> L_omega_scale;  // scale of LKJ prior for correlation of variance of m_0 (set to zero to ignore)
}

transformed data {
  real sigma_kappanu;
  vector[K] x_ss[M,L];          // sample sum of squares for each category (m) and group (l) in exposure
  for (m in 1:M)
    for (l in 1:L)
      x_ss[m,l] = (N[m,l] - 1) * x_sd[m,l];

  /* Scale for the prior of kappa/nu_0. In order to deal with input that does not contain observations
     (in which case n_each == 0), we set the minimum value for SD to 10. */
  sigma_kappanu = max(N) > 0 ? max(N) * 4 : 10;
}

parameters {
  // these are all shared across groups (same prior beliefs):
  vector<lower=K>[multiple_kappa_0 ? K : 1] kappa_0;        // prior pseudocount for category mu
  vector<lower=K + 1>[multiple_nu_0 ? K : 1] nu_0;          // prior pseudocount for category Sigma

  vector[K] m_0_param[m_0_known ? 0 : M];                 // prior mean of means
  vector<lower=0>[m_0_known ? 0 : K] m_0_tau;             // prior variances of m_0
  cholesky_factor_corr[m_0_known ? 0 : K] m_0_L_omega;    // prior correlations of variances of m_0 (in cholesky form)

  vector<lower=0>[K] S_0_param[S_0_known ? 0 : M];        // prior expected standard deviation

  real<lower=0, upper=1> lapse_rate;
}

transformed parameters {
  vector[K] m_0[M];                    // prior mean of means m_0
  vector[K] S_0[M];                    // prior expected standard deviation S_0

  // updated beliefs depend on input and group
  vector<lower=K>[multiple_kappa_0 ? K : 1] kappa_n[M,L];    // updated mean pseudocount
  vector<lower=K>[multiple_nu_0 ? K : 1] nu_n[M,L];          // updated sd pseudocount
  vector[K] m_n[M,L];                  // updated expected mean
  vector[K] S_n[M,L];                  // updated expected scatter matrix
  vector[K] t_scale[M,L];              // scale matrix of predictive t distribution

  simplex[M] p_test_conj[N_test];
  vector[M] log_p_test_conj[N_test];

  if (m_0_known) {
    m_0 = m_0_data;
  } else {
    m_0 = m_0_param;
  }
  if (S_0_known) {
    S_0 = S_0_data;
  } else {
    S_0 = S_0_param;
  }

  // update NIX parameters according to conjugate updating rules are taken from
  // Murphy (2007, p. XXX)
  for (cat in 1:M) {
    for (group in 1:L) {
      if (N[cat,group] > 0 ) {
        kappa_n[cat,group] = kappa_0 + N[cat,group];
        nu_n[cat,group] = nu_0 + N[cat,group];
        if (multiple_kappa_0) {
          m_n[cat,group] = (kappa_0 .* m_0[cat] + N[cat,group] * x_mean[cat,group]) ./ kappa_n[cat,group];
        } else {
          m_n[cat,group] = (kappa_0[1] * m_0[cat] + N[cat,group] * x_mean[cat,group]) / kappa_n[cat,group,1];
        }

        if (multiple_nu_0 && multiple_kappa_0) {
          S_n[cat,group] = sqrt(
            (nu_0 .* (S_0[cat] .* S_0[cat]) +
            x_ss[cat,group] +
            ((kappa_0 * N[cat,group]) ./ kappa_n[cat,group]) .*
            (m_0[cat] - x_mean[cat,group]) .* (m_0[cat] - x_mean[cat,group])) ./
            nu_n[cat,group]);
        } else {
          // NEED TO HANDLE CASES WHEN ONLY ONE OF MULTIPLE_NU_0 OR MULTIPLE_KAPPA_0 IS TRUE
          S_n[cat,group] = sqrt(
            (nu_0[1] * (S_0[cat] .* S_0[cat]) +
            x_ss[cat,group] +
            ((kappa_0[1] * N[cat,group]) / kappa_n[cat,group,1]) *
            (m_0[cat] - x_mean[cat,group]) .* (m_0[cat] - x_mean[cat,group])) /
            nu_n[cat,group, 1]);
        }


      } else {
        kappa_n[cat,group] = kappa_0;
        nu_n[cat,group] = nu_0;
        m_n[cat,group] = m_0[cat];
        S_n[cat,group] = S_0[cat];
      }

      t_scale[cat,group] = S_n[cat,group] .* (kappa_n[cat,group] + 1) ./
                                              (kappa_n[cat,group] .* (nu_n[cat,group] - K + 1));
    }
  }

  // compute category probabilities for each of the test stimuli
  for (j in 1:N_test) {
    int group;
    group = y_test[j];
    // calculate un-normalized log prob for each category
    for (cat in 1:M) {
      log_p_test_conj[j,cat] = multi_student_t_lpdf(x_test[j] |
                                              nu_n[cat,group] - K + 1,
                                              m_n[cat,group],
                                              t_scale[cat,group]);
    }
    // normalize and store actual probs in simplex
    p_test_conj[j] = exp(log_p_test_conj[j] - log_sum_exp(log_p_test_conj[j]));
  }
}

model {
  vector[M] lapsing_probs;

  lapsing_probs = rep_vector(lapse_rate / M, M);

  kappa_0 ~ normal(0, sigma_kappanu);
  nu_0 ~ normal(0, sigma_kappanu);

  /* Specifying prior for m_0:
     - If no scale for variances (tau) of m_0 is user-specified use weakly regularizing
       scale (5) for variances of mean.
     - If no scale for LKJ prior over correlation matrix of m_0 is user-specified use
       scale 1 to set uniform prior over correlation matrices. */
  if (!m_0_known) {
    m_0_tau ~ cauchy(0, tau_scale > 0 ? tau_scale : 5);
    m_0_L_omega ~ lkj_corr_cholesky(L_omega_scale > 0 ? L_omega_scale : 1);
    m_0_param ~ multi_normal_cholesky(rep_vector(0, K), diag_pre_multiply(m_0_tau, m_0_L_omega));
  }

  /* Specifying prior for components of S_0: */
  if (!S_0_known) {
      for (cat in 1:M) {
        tau_0_param[cat] ~ cauchy(0, tau_scale > 0 ? tau_scale : 10);
        L_omega_0_param[cat] ~ lkj_corr_cholesky(L_omega_scale > 0 ? L_omega_scale : 1);
      }
  }

  for (i in 1:N_test) {
    z_test_counts[i] ~ multinomial(p_test_conj[i] * (1-lapse_rate) + lapsing_probs);
  }

}

generated quantities {
  if (!m_0_known) {
    matrix[K,K] m_0_cor;
    matrix[K,K] m_0_cov;

    m_0_cor = multiply_lower_tri_self_transpose(m_0_L_omega);
    m_0_cov = quad_form_diag(m_0_cor, m_0_tau);
  }
}
