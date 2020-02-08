/*
 * Fit multinomial response data using a belief-updating model to infer prior
 * parameters (prior means and covariances can be inferred or provided by the
 * user). A normal-inverse-Wishart prior is used, and it's assumed
 * that the subject knows the true category labels of all stimuli during
 * exposure.
 *
 * This version has a lapse rate parameter (probability of random guessing).
 * The parameter can be inferred or provided by the user.
 *
 * Input is in the form of raw data points for observations
 *
 * Dave Kleinschmidt
 * Apr 2015
 *
 * Modified by
 * Florian Jaeger
 * June 2019
 *
 * Variable names follow Murphy (2012, p. 134):
 *
 * Category parameters: mu, sigma.
 * The mean (mu) and covariance matrix (sigma) of each category, indexed by beliefs.
 * E.g, mu_0 and sigma_0 are the mean and covariance under the prior beliefs, and mu_n
 * and sigma_n are the mean and covariance after n additional observations.
 *
 * Belief parameters: m, kappa, nu, S
 * The NIW describes the joint distribution of mu and sigma, as the product of a Normal
 * distribution over mu (expressing the uncertainty about the true category mean, mu)
 * multiplied with an Inverse-Wishart distribution over sigma (expressing the uncertainty
 * about the category covariance, sigma). The two components are linked through the sigma:
 *
 *    p(mu, sigma) = N(mu | m, 1/kappa * sigma) x IW(sigma | S, nu)
 */

data {
  int M;                               // number of categories
  int L;                               // number of grouping levels (e.g. subjects)
  int K;                               // number of features/cues

  /* Efficiency considerations: declaring n as matrix rather than array of int is supposedly faster
     (https://mc-stan.org/docs/2_18/stan-users-guide/indexing-efficiency-section.html). It also means
     that standard arithmetic functions can be directly applied to that matrix (unlike to more-than-1-
     dimensional arrays). */

  // Exposure data
  matrix[M,L] N;                       // number of observations per category (m) and subject (l)
  vector[K] x_mean[M,L];               // means for each category (m) and subject (l)
  cov_matrix[K] x_ss[M,L];             // sum of uncentered squares matrix for each category (m) and subject (l)

  // Test data
  int N_test;                          // number of test trials
  vector[K] x_test[N_test];            // locations (feature/cue values) of test trials
  int y_test[N_test];                  // group label of test trials
  int z_test_counts[N_test,M];         // responses for test trials

  /* The user can *optionally* specify the expected prior category mean (expected mu_0) and/or expected prior
     category covariance matrix (expected sigma_0). */
  int<lower = 0, upper = 1> prior_mu_known;    // are the prior expected means known and provided by the user?
  int<lower = 0, upper = 1> prior_sigma_known; // are the prior expected covariance matrices known and provided by the user?
  int<lower = 0, upper = 1> lapse_known;       // is the lapse rate known and provided by the user?

  vector[prior_mu_known ? 0 : K] mu_0_data[M]; // prior expected means (set to zero to ignore)
  cov_matrix[prior_sigma_known ? 0 : K] sigma_0_data[M];  // prior expected covariance matrices
  real<lower=0, upper=1> lapse_rate_data[lapse_known ? 0 : 1];

  /* These declarations could be made optional (tau_scale is only required if either mu_0 or tau_0 is unknown;
     L_omega_scale is only required if mu_0 or L_omega_0 are unknown). But it doesn't seem worth the additional
     decrease in coding transparency given that the model code already has if-statements in place to do this
     work */
  real<lower=0> tau_scale;             // scale of Cauchy prior for variances of mu_0 (set to zero to ignore)
  real<lower=0> L_omega_scale;         // scale of LKJ prior for correlation of variance of mu_0 (set to zero to ignore)
}

transformed data {
  real sigma_kappanu;                  // variance of Normal prior that determines likely range of kappa_0 and nu_0

  sigma_kappanu = max(N) > 0 ? max(N) * 4 : 10;
}

parameters {
  /* These parameters are all shared across subjects and exposure conditions (i.e., everyone is assumed to
     have the same prior beliefs). */
  real<lower=K> kappa_0;               // prior pseudocount for mean
  real<lower=K> nu_0;                  // prior pseudocount for sd

  /* If mu_0 is known, set these parameters for the prior distribution of the mean to have zero dimensionality. */
  vector[!prior_mu_known ? 0 : K] m_0_param[M];                // prior expected means (set to zero to ignore)
  vector<lower=0>[!prior_mu_known ? 0 : K] m_0_tau;            // prior variances of mu_0
  cholesky_factor_corr[!prior_mu_known ? 0 : K] m_0_L_omega;   // prior correlations of variances of mu_0

  /* If sigma_0 is known, set parameters for the prior distribution of the standard deviations and correlation
     matrices to have zero dimensionality. */
  vector<lower=0>[!prior_sigma_known ? 0 : K] tau_0[M];          // prior standard deviations of scatter matrixs_0
  cholesky_factor_corr[!prior_sigma_known ? 0 : K] L_omega_0[M]; // prior correlations of variances of scatter matrix s_0

  real<lower=0, upper=1> lapse_rate_param[!lapse_known ? 0 : 1];
}

transformed parameters {
  /* Define the parameters for the prior distribution of means, standard deviations, and correlation matrices,
     depending on whether they were provided by the user (and are thus 'known') or not. */
  vector[K] m_0[M];                    // prior mean of Normal over category means
  cov_matrix[K] s_0[M];                // prior scatter matrix of IW
  real<lower=0, upper=1> lapse_rate;

  // Updated beliefs depend on input/subject
  real<lower=K> kappa_n[M,L];          // updated mean pseudocount
  real<lower=K> nu_n[M,L];             // updated sd pseudocount
  vector[K] m_n[M,L];                  // updated mean of Normal over category means
  cov_matrix[K] s_n[M,L];              // updated scatter matrix of IW
  cov_matrix[K] t_scale_n[M,L];        // updated scale matrix of predictive t distribution

  simplex[M] p_test_conj[N_test];
  vector[M] log_p_test_conj[N_test];

  if (prior_mu_known) {
    m_0 = mu_0_data;
  } else {
    m_0 = m_0_param;
  }

  if (lapse_known) {
    lapse_rate = lapse_rate_data[1];
  } else {
    lapse_rate = lapse_rate_param[1];
  }

  // update NIW parameters according to conjugate updating rules are taken from
  // Murphy (2007, p. 136)
  for (cat in 1:M) {
    if (prior_sigma_known) {
      // Convert the expected covariance matrix of each category into the scatter matrix s_0.
      s_0[cat] = sigma_0_data[cat] * (nu_0 - K - 1);
    } else {
      // Get sigma_0 from its components: correlation matrix and vector of standard deviations
      s_0[cat] = quad_form_diag(multiply_lower_tri_self_transpose(L_omega_0[cat]), tau_0[cat]);
    }

    for (subj in 1:L) {
      if (N[cat,subj] > 0 ) {
        kappa_n[cat,subj] = kappa_0 + N[cat,subj];
        nu_n[cat,subj] = nu_0 + N[cat,subj];
        m_n[cat,subj] = (kappa_0 * m_0[cat] + N[cat,subj] * x_mean[cat,subj]) /
                        kappa_n[cat,subj];
        s_n[cat,subj] = s_0[cat] +
                        x_ss[cat,subj] +
                        kappa_0 * m_0[cat] * m_0[cat]' -
                        kappa_n[cat,subj] * m_n[cat,subj] * m_n[cat,subj]';
      } else {
        kappa_n[cat,subj] = kappa_0;
        nu_n[cat,subj] = nu_0;
        m_n[cat,subj] = m_0[cat];
        s_n[cat,subj] = s_0[cat];
      }

      t_scale_n[cat,subj] = s_n[cat,subj] * (kappa_n[cat,subj] + 1) /
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
                                              m_n[cat,subj],
                                              t_scale_n[cat,subj]);
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

  /* If prior mu is NOT known, specifying prior for mu:
     - If no scale for variances (tau) of mu_0 is user-specified use weakly regularizing
       scale (5) for variances of mean.
     - If no scale for LKJ prior over correlation matrix of mu_0 is user-specified use
       scale 1 to set uniform prior over correlation matrices. */
  if (!prior_mu_known) {
    m_0_tau ~ cauchy(0, tau_scale > 0 ? tau_scale : 5);
    m_0_L_omega ~ lkj_corr_cholesky(L_omega_scale > 0 ? L_omega_scale : 1);
    m_0 ~ multi_normal_cholesky(rep_vector(0, K), diag_pre_multiply(m_0_tau, m_0_L_omega));
  }

  /* If prior sigma is NOT known, set priors for tau_0 and L_omega_0. */
  if (!prior_sigma_known) {
    for (cat in 1:M) {
      tau_0[cat] ~ cauchy(0, tau_scale > 0 ? tau_scale : 10);
      L_omega_0[cat] ~ lkj_corr_cholesky(L_omega_scale > 0 ? L_omega_scale : 1);
    }
  }

  for (i in 1:N_test) {
    z_test_counts[i] ~ multinomial(p_test_conj[i] * (1-lapse_rate) + lapsing_probs);
  }

}

generated quantities {
  matrix[K,K] m_0_cor;
  matrix[K,K] m_0_cov;

  m_0_cor = multiply_lower_tri_self_transpose(m_0_L_omega);
  m_0_cov = quad_form_diag(m_0_cor, ms_0_tau);
}
