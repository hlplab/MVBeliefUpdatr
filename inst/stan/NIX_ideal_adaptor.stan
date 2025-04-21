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
  int M;                        // number of categories
  int L;                        // number of groups (conditions or subjects)

  matrix[M,L] N_exposure;             // number of training observations
  matrix[M,L] x_mean_exposure;        // mean of training observations
  matrix[M,L] x_sd_exposure;          // sample standard deviation of training obs

  int N_test;                         // number of test trials
  array[N_test] real x_test;          // locations of test trials
  array[N_test] int y_test;           // group labels for test trials
  array[N_test,M] int z_test_counts;  // responses for test trials

  // NOT YET USED
  int<lower=0, upper=1> lapse_rate_known;
  int<lower=0, upper=1> mu_0_known;
  int<lower=0, upper=1> sigma_0_known;
}

transformed data {
  matrix[M,L] x_ss_exposure = (N_exposure - 1) .* x_sd_exposure;

  real<lower=0> sigma_kappanu = max(to_array_1d(N_exposure .* 4)); // scale for the prior of kappa/nu_0 (at least 10)
}

parameters {
  // these are all shared across groups (same prior beliefs):
  real<lower=0> kappa_0;       // prior pseudocount for mean
  real<lower=0> nu_0;          // prior pseudocount for sd
  array[M] real m_0;                 // prior expected mean
  array[M] real<lower=0> S_0;        // prior expected standard deviation
  real<lower=0, upper=1> lapse_rate;
}

transformed parameters {
  // updated beliefs depend on input/group
  array[M,L] real m_n;                 // updated expected mean
  array[M,L] real<lower=0> kappa_n;    // updated mean pseudocount
  array[M,L] real<lower=0> S_n;        // updated expected sd
  array[M,L] real<lower=0> nu_n;       // updated sd pseudocount
  array[M,L] real<lower=0> t_scale;    // scale parameter of predictive t distribution

  /* Assuming unifom bias, so that lapsing_prob = probability of each category prior to
     taking into account stimulus is 1/M
  */
  vector[M] lapsing_probs = rep_vector(lapse_rate / M, M);

  array[N_test] simplex[M] p_test_conj;
  array[N_test] vector[M] log_p_test_conj;

  // update NIX2 parameters according to conjuate updating rules are taken from
  // Murphy (2007)
  for (cat in 1:M) {
    for (group in 1:L) {
      kappa_n[cat,group] = kappa_0 + N_exposure[cat,group];
      nu_n[cat,group] = nu_0 + N_exposure[cat,group];
      m_n[cat,group] = (m_0[cat] * kappa_0 + x_mean_exposure[cat,group] * N_exposure[cat,group]) / kappa_n[cat,group];
      S_n[cat,group] = sqrt((nu_0*S_0[cat]^2 +
                                 x_ss_exposure[cat,group] +
                                 (N_exposure[cat,group] * kappa_0) / (kappa_n[cat,group]) *
                                   (m_0[cat] - x_mean_exposure[cat,group])^2
                                 ) /
                                nu_n[cat,group]);
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

  m_0 ~ normal(0, 100);
  S_0 ~ uniform(0, 100);

  for (i in 1:N_test) {
    z_test_counts[i] ~ multinomial(p_test_conj[i] * (1-lapse_rate) + lapsing_probs);
  }

}

generated quantities {

}
