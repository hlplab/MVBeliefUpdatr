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
    int L;                        // number of exposure groups (e.g. subjects)
    int K;                        // number of features

    /* Efficiency considerations: declaring n as matrix rather than array of int is supposedly faster
    (https://mc-stan.org/docs/2_18/stan-users-guide/indexing-efficiency-section.html). It also means
    that standard arithmetic functions can be directly applied to that matrix (unlike to more-than-1-
                                                                                 dimensional arrays). */
    matrix[M,L] N;                // number of observations per category (M) and exposure group (L)
    vector[K] x_mean[M,L];        // means for each category (M) and exposure group (L)
    cov_matrix[K] x_ss[M,L];      // sum of uncentered squares matrix for each category (M) and group (L)

    int N_test;                   // number of unique combinations of test locations & exposure groups
    vector[K] x_test[N_test];     // locations (in cue space) of test trials
    int y_test[N_test];           // exposure group indicator of test trials
    int z_test_counts[N_test,M];  // responses for test trials (for each category M)

    int<lower=0, upper=1> lapse_rate_known;
    int<lower=0, upper=1> mu_0_known;
    int<lower=0, upper=1> Sigma_0_known;
    real lapse_rate_data[lapse_rate_known ? 1 : 0];                         // optional: user provided lapse_rate
    vector[mu_0_known ? K : 0] mu_0_data[mu_0_known ? M : 0];               // optional: user provided expected mu_0 (prior category means)
    cov_matrix[Sigma_0_known ? K : 0] Sigma_0_data[Sigma_0_known ? M : 0];  // optional: user provided expected Sigma_0 (prior category covariance matrices)

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

  /* Scale for the prior of kappa/nu_0. In order to deal with input that does not contain observations
  (in which case n_each == 0), we set the minimum value for SD to 10. */
    sigma_kappanu = max(N) > 0 ? max(N) * 4 : 10;
}

parameters {
  // these are all shared across groups (same prior beliefs):
    real<lower=K> kappa_0;                  // prior pseudocount for category mu
    real<lower=K + 1> nu_0;                 // prior pseudocount for category Sigma

    real<lower=0, upper=1> lapse_rate_param[lapse_rate_known ? 0 : 1];

    vector[K] m_0_param[mu_0_known ? 0 : M];                // prior mean of means
    vector<lower=0>[mu_0_known ? 0 : K] m_0_tau;            // prior variances of m_0
    cholesky_factor_corr[mu_0_known ? 0 : K] m_0_L_omega;   // prior correlations of variances of m_0 (in cholesky form)

    vector<lower=0>[K] tau_0_param[Sigma_0_known ? 0 : M];          // standard deviations of prior scatter matrix S_0
    cholesky_factor_corr[K] L_omega_0_param[Sigma_0_known ? 0 : M]; // correlation matrix of prior scatter matrix S_0 (in cholesky form)
}

transformed parameters {
  real lapse_rate;
  vector[K] m_0[M];                    // prior mean of means m_0
  cov_matrix[K] S_0[M];                // prior scatter matrix S_0

  // updated beliefs depend on input and group
  real<lower=K> kappa_n[M,L];          // updated mean pseudocount
  real<lower=K> nu_n[M,L];             // updated sd pseudocount
  vector[K] m_n[M,L];                  // updated expected mean
  cov_matrix[K] S_n[M,L];              // updated expected scatter matrix
  cov_matrix[K] t_scale[M,L];          // scale matrix of predictive t distribution

  simplex[M] p_test_conj[N_test];
  vector[M] log_p_test_conj[N_test];

  if (lapse_rate_known) {
    lapse_rate = lapse_rate_data[1];
  } else {
    lapse_rate = lapse_rate_param[1];
  }
  if (mu_0_known) {
    m_0 = mu_0_data;
  } else {
    m_0 = m_0_param;
  }
  if (Sigma_0_known) {
   // Get S_0 from expected Sigma given nu_0
   for (cat in 1:M) {
      S_0[cat] = Sigma_0_data[cat] * (nu_0 - K - 1);
   }
  }

  // update NIW parameters according to conjugate updating rules are taken from
  // Murphy (2007, p. 136)
  for (cat in 1:M) {
    if (!Sigma_0_known) {
      // Get S_0 from its components: correlation matrix and vector of standard deviations
      S_0[cat] = quad_form_diag(multiply_lower_tri_self_transpose(L_omega_0_param[cat]), tau_0_param[cat]);
    }
    for (group in 1:L) {
      if (N[cat,group] > 0 ) {
        kappa_n[cat,group] = kappa_0 + N[cat,group];
        nu_n[cat,group] = nu_0 + N[cat,group];
        m_n[cat,group] = (kappa_0 * m_0[cat] + N[cat,group] * x_mean[cat,group]) /
          kappa_n[cat,group];
        S_n[cat,group] = S_0[cat] +
          x_ss[cat,group] +
          kappa_0 * m_0[cat] * m_0[cat]' -
                        kappa_n[cat,group] * m_n[cat,group] * m_n[cat,group]';
      } else {
        kappa_n[cat,group] = kappa_0;
        nu_n[cat,group] = nu_0;
        m_n[cat,group] = m_0[cat];
        S_n[cat,group] = S_0[cat];
      }

      t_scale[cat,group] = S_n[cat,group] * (kappa_n[cat,group] + 1) /
        (kappa_n[cat,group] * (nu_n[cat,group] - K + 1));
    }
  }

  // compute posterior probabilities of all categories for each of the test stimuli
  for (j in 1:N_test) {
    int group;
    group = y_test[j];
    /* calculate un-normalized log posterior probability for each category. Under the assumption
       of uniform prior probabilities for each category, the log probabilities identical to the
       normalized log likelihoods. If we ever were to change this assumption, we'd have to add
       the log prior probabilities of categories here. */
    for (cat in 1:M) {
      // log likelihood of the test stimulus given the category
      log_p_test_conj[j,cat] = multi_student_t_lpdf(x_test[j] |
                                                    nu_n[cat,group] - K + 1,
                                                    m_n[cat,group],
                                                    t_scale[cat,group]);
    }
    // normalize to get actual posterior probs in simplex
    p_test_conj[j] = exp(log_p_test_conj[j] - log_sum_exp(log_p_test_conj[j]));
  }
}

model {
  vector[M] lapsing_probs;

  /* Assuming unifom bias, so that lapsing_prob = probability of each category prior to
     taking into account stimulus is 1/M */
  lapsing_probs = rep_vector(lapse_rate / M, M);

  kappa_0 ~ normal(0, sigma_kappanu);
  nu_0 ~ normal(0, sigma_kappanu);

  /* Specifying prior for m_0:
     - If no scale for variances (tau) of m_0 is user-specified use weakly regularizing
       scale (5) for variances of mean.
     - If no scale for LKJ prior over correlation matrix of m_0 is user-specified use
       scale 1 to set uniform prior over correlation matrices. */
  if (!mu_0_known) {
      m_0_tau ~ cauchy(0, tau_scale > 0 ? tau_scale : 5);
      m_0_L_omega ~ lkj_corr_cholesky(L_omega_scale > 0 ? L_omega_scale : 1);
      m_0_param ~ multi_normal_cholesky(rep_vector(0, K), diag_pre_multiply(m_0_tau, m_0_L_omega));
  }

  // Specifying prior for components of S_0:
  if (!Sigma_0_known) {
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
  /* Compute and store pointwise log-likelihood, in order to allow computation of LOOIC.
     Doing so in generated quantities block, following help(rstan::loo). Note that currently,
     each unique combination of test location and exposure group is treated as an observation
     (rather than each individual response). In the future, this should probably changed to
     treated each individual response as an observation. */
  vector[M] lapsing_probs = rep_vector(lapse_rate / M, M);
  vector[N_test] log_lik;
  real log_lik_sum = 0;
  for (n in 1:N_test) {
    log_lik[n] = multinomial_lpmf(z_test_counts[n] | p_test_conj[n] * (1-lapse_rate) + lapsing_probs);
    log_lik_sum += log_lik[n];
  }

  if (!mu_0_known) {
    matrix[K,K] m_0_cor;
    matrix[K,K] m_0_cov;

    m_0_cor = multiply_lower_tri_self_transpose(m_0_L_omega);
    m_0_cov = quad_form_diag(m_0_cor, m_0_tau);
  }
}
