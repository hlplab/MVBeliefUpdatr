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
* May 2025
*/

data {
  int M;                                                   // number of categories
  int L;                                                   // number of exposure groups (e.g. subjects)
  int K;                                                   // number of features

  array[M,L] int<lower=0> N;                               // number of observations per category (M) and exposure group (L)
  array[M,L] vector[K] x_mean;                             // means for each category (M) and exposure group (L)
  array[M,L] cov_matrix[K] x_ss;                           // sum of uncentered squares matrix for each category (M) and group (L)

  int N_test;                                              // number of unique combinations of test locations & exposure groups
  array[N_test] vector[K] x_test;                          // locations (in cue space) of test trials
  array[N_test] int y_test;                                // exposure group indicator of test trials
  array[N_test,M] int z_test_counts;                       // responses for test trials (for each category M)

  int<lower=0, upper=1> lapse_rate_known;
  int<lower=0, upper=1> mu_0_known;
  int<lower=0, upper=1> Sigma_0_known;
  array[lapse_rate_known ? 1 : 0] real<lower=0, upper=1> lapse_rate_data;       // optional: user provided lapse_rate
  array[mu_0_known ? M : 0] vector[mu_0_known ? K : 0] mu_0_data;               // optional: user provided expected mu_0 (prior category means)
  array[Sigma_0_known ? M : 0] cov_matrix[Sigma_0_known ? K : 0] Sigma_0_data;  // optional: user provided expected Sigma_0 (prior category covariance matrices)

  /* For now, this script assumes that the observations (cue vectors) are centered. The prior
     mean of m_0 is set to 0. In the future, one could automatically adjust the location based
     on the input data (e.g., by placing the mean of the Normal prior in the center of the
     exposure data).
  */
  vector<lower=0>[mu_0_known ? 0 : K] tau_scale;           // separate taus for each of the K features to capture that features can be on separate scales
  real<lower=0> L_omega_eta;                                // eta of LKJ prior for correlation of variance of m_0

  int<lower=0, upper=1> split_loglik_per_observation;
}

transformed data {
  /* Scale for the prior of kappa/nu_0. In order to deal with input that does not contain observations
     (in which case n_each == 0), we set the minimum value for SD to 10.
  */
  real<lower=0> sigma_kappanu = min(max(to_array_1d(N)) * 4, 10);
  vector[K] m_0_mu = rep_vector(0, K);         // center of prior of m_0
}

parameters {
  // these are all shared across groups (same prior beliefs):
  real<lower=K> kappa_0;                                   // prior pseudocount for category mu
  real<lower=K+1> nu_0;                                    // prior pseudocount for category Sigma

  array[lapse_rate_known ? 0 : 1] real<lower=0, upper=1> lapse_rate_param;

  array[mu_0_known ? 0 : M] vector[K] m_0_param;           // prior mean of means
  vector<lower=0>[mu_0_known ? 0 : K] m_0_tau;             // prior variances of m_0
  cholesky_factor_corr[mu_0_known ? 0 : K] m_0_L_omega;    // prior correlations of variances of m_0 (in cholesky form)

  array[Sigma_0_known ? 0 : M] vector<lower=0>[K] tau_0_param;          // standard deviations of prior scatter matrix S_0
  array[Sigma_0_known ? 0 : M] cholesky_factor_corr[K] L_omega_0_param; // correlation matrix of prior scatter matrix S_0 (in cholesky form)
}

transformed parameters {
  array[M] vector[K] m_0 = mu_0_known ? mu_0_data : m_0_param;                    // prior mean of means m_0
  array[M] cholesky_factor_cov[K] L_S_0;                                          // Cholesky factor of prior scatter matrix S_0
  real lapse_rate = lapse_rate_known ? lapse_rate_data[1] : lapse_rate_param[1];  // lapse rate
  // Assuming unifom bias, so that lapsing_prob = probability of each category prior to
  // taking into account stimulus is 1/M
  vector[M] lapsing_probs = rep_vector(lapse_rate / M, M);

  // updated beliefs depend on input and group
  array[M,L] real<lower=K> kappa_n;             // updated mean pseudocount
  array[M,L] real<lower=K> nu_n;                // updated sd pseudocount
  array[M,L] vector[K] m_n;                     // updated expected mean
  array[M,L] cholesky_factor_cov[K] L_S_n;      // Cholesky factor of updated expected scatter matrix
  array[M,L] cholesky_factor_cov[K] L_t_scale;  // Cholesky factor of scale matrix of predictive t distribution

  array[N_test] simplex[M] p_test_conj;
  array[N_test] vector[M] log_p_test_conj;

  // update NIW parameters according to conjugate updating rules are taken from
  // Murphy (2007, p. 136)
  for (cat in 1:M) {
    /* Get S_0 from expected Sigma given nu_0

       Using Cholesky factor representation of S_0 in order to avoid issues with numerical
       instability that would otherwise arise for S_0s with large values (leading to exceptions
       caused by non-symmetrical S_0). The use of Cholesky factors should also make the
       computation faster.

       The idea was provided by Bob Carpenter (thanks!) and implemented below with the help
       of ChatGPT.
       https://discourse.mc-stan.org/t/exception-s-n-cov-matrix-not-symmetric-for-large-value-elements/39184/2
    */
    L_S_0[cat] = Sigma_0_known ? cholesky_decompose(Sigma_0_data[cat] * (nu_0 - K - 1)) :
                            diag_pre_multiply(tau_0_param[cat], L_omega_0_param[cat]);
    for (group in 1:L) {
      if (N[cat,group] > 0 ) {
        kappa_n[cat,group] = kappa_0 + N[cat,group];
        nu_n[cat,group] = nu_0 + N[cat,group];
        m_n[cat,group] = (kappa_0 * m_0[cat] + N[cat,group] * x_mean[cat,group]) /
          kappa_n[cat,group];
        L_S_n[cat, group] = cholesky_decompose(
          symmetrize_from_lower_tri(
            multiply_lower_tri_self_transpose(L_S_0[cat]) +
            x_ss[cat, group] +
            kappa_0 * m_0[cat] * m_0[cat]' -
            kappa_n[cat, group] * m_n[cat, group] * m_n[cat, group]'
          )
        );

      } else {
        kappa_n[cat,group] = kappa_0;
        nu_n[cat,group] = nu_0;
        m_n[cat,group] = m_0[cat];
        L_S_n[cat,group] = L_S_0[cat];
      }

      // Instead of computing t_scale as in Murphy, we calculate the Cholesky factor of that
      // scale, and then use the Cholesky-based form of the multivariate T-density below:
      L_t_scale[cat, group] = L_S_n[cat, group] * sqrt((kappa_n[cat, group] + 1) /
                               (kappa_n[cat, group] * (nu_n[cat, group] - K + 1)));
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
      // Use Cholesky form of multivariate T-density
      log_p_test_conj[j,cat] =
          multi_student_t_cholesky_lpdf(
              x_test[j] |
              nu_n[cat,group] - K + 1,
              m_n[cat,group],
              L_t_scale[cat,group]);
    }
    // normalize to get actual posterior probs in simplex
    p_test_conj[j] = exp(log_p_test_conj[j] - log_sum_exp(log_p_test_conj[j]));
  }
}

model {
  kappa_0 ~ normal(0, sigma_kappanu);
  nu_0 ~ normal(0, sigma_kappanu);

  if (!mu_0_known) {
      m_0_tau ~ cauchy(0, tau_scale);
      m_0_L_omega ~ lkj_corr_cholesky(L_omega_eta);
      m_0_param ~ multi_normal_cholesky(m_0_mu, diag_pre_multiply(m_0_tau, m_0_L_omega));
  }

  // Specifying prior for components of S_0:
  if (!Sigma_0_known) {
    for (cat in 1:M) {
      tau_0_param[cat] ~ cauchy(0, tau_scale);
      L_omega_0_param[cat] ~ lkj_corr_cholesky(L_omega_eta);
    }
  }

  for (i in 1:N_test) {
    z_test_counts[i] ~ multinomial(p_test_conj[i] * (1 - lapse_rate) + lapsing_probs);
  }
}

generated quantities {
  /* Compute and store pointwise log-likelihood, in order to allow computation of LOOIC.
     Doing so in generated quantities block, following help(rstan::loo). Note that currently,
     each unique combination of test location and exposure group is treated as an observation
     (rather than each individual response).
  */
  real log_lik_sum = 0;
  vector[N_test] log_lik;

  for (n in 1:N_test) {
    log_lik[n] = multinomial_lpmf(z_test_counts[n] | p_test_conj[n] * (1-lapse_rate) + lapsing_probs);
    log_lik_sum += log_lik[n];
  }

  /* If requested by user:
     Compute and store pointwise log-likelihood, in order to allow computation of LOOIC. This is
     done in the generated quantities block, following help(rstan::loo). Note that each unqiue
     combination of test location and exposure group is treated as an observation in z_test_counts.
     So calculating pointwise log-likelihoods based on the vectorized z_test_counts, which would be
     most straightforward, has several downsides:

     1) it makes the pointwise log-likelihoods, and differences between the individual pointwise
        log-likelihoods, hard to interpret (e.g., because test locations and/or groups differ in
        the number of observations, which affects the log-likelihoods).

      2) it makes it impossible to compare across different models that are fit on the same under-
         lying data if they differ in the number of unique combinations of test locations and
         groups (and thus the number of rows of z_test_counts). This can happen, for instance, when
         one compares models over different transformation of the cues, which can result in different
         numbers of unique test locations.

     Here, log-likelihoods are thus calculated for each individual response. This requires reformatting
     the information z_test_counts.

     The downside of this revised approach is that it makes the object very large to store, and slow to
     work with.
  */
  vector[split_loglik_per_observation ? sum(to_array_1d(z_test_counts)) : 0] log_lik_split;
  if (split_loglik_per_observation == 1) {
    int idx = 1;

    for (n in 1:N_test) {
      for (m in 1:M) {
        for (i in 1:z_test_counts[n, m]) {
          log_lik_split[idx] = log(p_test_conj[n, m] * (1 - lapse_rate) + lapsing_probs[m]);
          idx += 1;
        }
      }
    }
  }

  matrix[mu_0_known ? 0 : K,mu_0_known ? 0 : K] m_0_cor;
  matrix[mu_0_known ? 0 : K,mu_0_known ? 0 : K] m_0_cov;
  if (!mu_0_known) {
    m_0_cor = multiply_lower_tri_self_transpose(m_0_L_omega);
    m_0_cov = quad_form_diag(m_0_cor, m_0_tau);
  }

  // Get S_0 and S_n from their Cholesky factors
  array[M] cov_matrix[K] S_0;         // prior scatter matrix S_0
  array[M,L] cov_matrix[K] S_n;              // Updated expected scatter matrix
  for (cat in 1:M) {
    S_0[cat] = multiply_lower_tri_self_transpose(L_S_0[cat]);
    for (group in 1:L) {
      S_n[cat,group] = multiply_lower_tri_self_transpose(L_S_n[cat,group]);
    }
  }
}
