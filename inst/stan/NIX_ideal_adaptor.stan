/*
 * Fit multinomial response data using a belief-updating model to infer prior
 * parameters.  A normal-inverse-chi-squared prior is used, and it's assumed
 * that the subject knows the true labels of all the input stimuli (e.g.,
 * doesn't model any uncertainty in classification.
 *
 * This version has a lapse rate parameter (probability of random guessing)
 *
 * Input is in the form of the sufficient statistics of each category during
 * exposure
 *
 * Dave Kleinschmidt
 * Feb 2017
 */

data {
  int M;                              // number of categories
  int L;                              // number of groups (conditions or subjects)

  matrix[M,L] N_exposure;             // number of training observations
  matrix[M,L] x_mean_exposure;        // mean of training observations
  matrix[M,L] x_sd_exposure;          // sample standard deviation of training obs

  int N_test;                         // number of test trials
  array[N_test] real x_test;          // locations of test trials
  array[N_test] int y_test;           // group labels for test trials
  array[N_test,M] int z_test_counts;  // responses for test trials

  int<lower=0, upper=1> lapse_rate_known;
  int<lower=0, upper=1> mu_0_known;
  int<lower=0, upper=1> Sigma_0_known;
  array[lapse_rate_known ? 1 : 0] real<lower=0, upper=1> lapse_rate_data;   // optional: user provided lapse_rate
  array[mu_0_known ? M : 0] real mu_0_data;                                 // optional: user provided expected mu_0 (prior category means) in space of affine transformation defined by INV_SCALE^-1 and shift
  array[Sigma_0_known ? M : 0] real Sigma_0_data;                           // optional: user provided expected Sigma_0 (prior category variances) in space of affine transformation defined by INV_SCALE^-1 and shift

  /* For now, this script assumes that the observations (cue vectors) are centered. The prior
     mean of m_0 is set to 0. In the future, one could automatically adjust the location based
     on the input data (e.g., by placing the mean of the Normal prior in the center of the
     exposure data).
  */
  real<lower=0> tau_scale;            // scale for the prior over means and variances

  /* The data above are assumed to have been transformed by a sensible affine transformation f(x) = A(x - center)
     (e.g., by scaling or whitening the data) in order to improve numerical stability. A_inv and center are used
     below in order to transform the parameters of interest back into the original space.
  */
  real INV_SCALE;                     // Inverse of transformation (must be invertible) that was applied to data
  real shift;                         // Center of affine transformation that was applied to data

  int<lower=0, upper=1> split_loglik_per_observation;
}

transformed data {
  matrix[M,L] x_ss_exposure = (N_exposure - 1) .* x_sd_exposure;

  real<lower=0> sigma_kappanu = max(to_array_1d(N_exposure .* 4)); // scale for the prior of kappa/nu_0 (at least 10)
  real m_0_mu = 0;                                                 // center of prior of m_0
}

parameters {
  /* The strength of the prior beliefs is currently assumed to be shared across groups and categories.
     Any of these assumptions could be revisited in the future. For example, different groups could differ in
     the strength of their prior beliefs because of differences in prior experience. Different categories could
     be associated with weaker or stronger prior beliefs because of differences in the amount of exposure
     (overall frequency of the category) or because differences in e.g., the cross-talker variability in the
     realization of the category.
  */
  real<lower=0> kappa_0;             // prior pseudocount for mean
  real<lower=2> nu_0;                // prior pseudocount for sd

  array[lapse_rate_known ? 0 : 1] real<lower=0, upper=1> lapse_rate_param;
  array[mu_0_known ? 0 : M] real m_0_param;     // prior expected means
  array[mu_0_known ? 0 : M] real m_0_tau;       // prior variances of expected means (m_0)
  array[Sigma_0_known ? 0 : M] real S_0_param;  // prior
}

transformed parameters {
  array[M] real m_0 = mu_0_known ? mu_0_data : m_0_param;                         // prior expected means m_0
  array[M] real S_0;                                                              // prior scatter matrix S_0
  real lapse_rate = lapse_rate_known ? lapse_rate_data[1] : lapse_rate_param[1];  // lapse rate

  /* Assuming unifom bias, so that lapsing_prob = probability of each category prior to
     taking into account stimulus is 1/M
  */
  vector[M] lapsing_probs = rep_vector(lapse_rate / M, M);

  // updated beliefs depend on input/group
  array[M,L] real<lower=0> kappa_n;    // updated mean pseudocount
  array[M,L] real<lower=0> nu_n;       // updated sd pseudocount
  array[M,L] real m_n;                 // updated expected mean
  array[M,L] real<lower=0> S_n;        // updated expected sd
  array[M,L] real<lower=0> t_scale;    // scale parameter of predictive t distribution

  array[N_test] simplex[M] p_test_conj;
  array[N_test] vector[M] log_p_test_conj;

  // update NIX2 parameters according to conjuate updating rules are taken from
  // Murphy (2007, p. 136). NOTE: Murphy reports E[sigma^2 | D] = nu / (nu - 2) * S,
  // but both comparison to the NIW and Google searches suggest that the * nu factor is an error (and so it is omitted here)
  // alternative, we would use S_0[cat] = Sigma_0_known ? Sigma_0_data[cat] * (nu_0 - 2) / nu_0: S_0_param[cat];
  for (cat in 1:M) {
    S_0[cat] = Sigma_0_known ? Sigma_0_data[cat] * (nu_0 - 2) : S_0_param[cat];
    for (group in 1:L) {
      if (N_exposure[cat,group] > 0 ) {
      kappa_n[cat,group] = kappa_0 + N_exposure[cat,group];
      nu_n[cat,group] = nu_0 + N_exposure[cat,group];
      m_n[cat,group] = (m_0[cat] * kappa_0 + x_mean_exposure[cat,group] * N_exposure[cat,group]) / kappa_n[cat,group];
      S_n[cat,group] = sqrt((nu_0*S_0[cat]^2 +
                                 x_ss_exposure[cat,group] +
                                 (N_exposure[cat,group] * kappa_0) / (kappa_n[cat,group]) *
                                   (m_0[cat] - x_mean_exposure[cat,group])^2
                                 ) /
                                nu_n[cat,group]);
      } else {
        kappa_n[cat,group] = kappa_0;
        nu_n[cat,group] = nu_0;
        m_n[cat,group] = m_0[cat];
        S_n[cat,group] = S_0[cat];
      }

      t_scale[cat,group] = S_n[cat,group] * sqrt((kappa_n[cat,group] + 1) / kappa_n[cat,group]);
    }
  }

  // compute category probabilities for each of the test stimuli
  for (j in 1:N_test) {
    int group;
    group = y_test[j];
    // calculate un-normalized log prob for each category
    for (cat in 1:M) {
      log_p_test_conj[j,cat] = student_t_lpdf(x_test[j] |
                                              nu_n[cat,group],
                                              m_n[cat,group],
                                              t_scale[cat,group]);
    }
    // normalize and store actual probs in simplex
    p_test_conj[j] = exp(log_p_test_conj[j] - log_sum_exp(log_p_test_conj[j]));
  }
}

model {
  kappa_0 ~ normal(0, sigma_kappanu);
  nu_0 ~ normal(0, sigma_kappanu);

  if (!mu_0_known) {
      m_0_tau ~ cauchy(0, tau_scale);
      m_0_param ~ normal(m_0_mu, m_0_tau);
  }

  if (!Sigma_0_known) {
      S_0_param ~ cauchy(0, tau_scale);
  }

  for (n in 1:N_test) {
    z_test_counts[n] ~ multinomial(p_test_conj[n] * (1-lapse_rate) + lapsing_probs);
  }

}

generated quantities {
  real log_lik_sum = 0;
  vector[N_test] log_lik;

  for (n in 1:N_test) {
    log_lik[n] = multinomial_lpmf(z_test_counts[n] | p_test_conj[n] * (1-lapse_rate) + lapsing_probs);
    log_lik_sum += log_lik[n];
  }

  // If requested by user:
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

  // Inverse transform parameters back into the original space
  array[M] real m_0_original;
  array[M] real S_0_original;
  array[M,L] real m_n_original;
  array[M,L] real S_n_original;
  for (cat in 1:M) {
    /* The validity of back-transforming a NIW by back-transforming its m and S parameters (kappa and nu
       remain unchanged) was confirmed via ChatGPT (https://chatgpt.com/c/683b3ed6-4924-800c-8f22-e61680f3360f)
    */
    // Since we define the affine transform as f(x) = SCALE(x + shift), rather than f(x) = SCALE * x + shift:
    m_0_original[cat] = INV_SCALE * m_0[cat] - shift;
    // Since the inverse-Wishart distribution is affine-invariant under congruence transformations:
    S_0_original[cat] = INV_SCALE^2 * S_0[cat];
    for (group in 1:L) {
      m_n_original[cat,group] = INV_SCALE * m_n[cat,group] - shift;
      S_n_original[cat,group] = INV_SCALE^2 * S_n[cat,group];
    }
  }
}
