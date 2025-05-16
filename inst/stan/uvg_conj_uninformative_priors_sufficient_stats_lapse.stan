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
  int m;                        // number of categories
  int l;                        // number of groups (conditions or subjects)

  matrix[m,l] n;                // number of training observations
  matrix[m,l] xbar;             // mean of training observations
  matrix[m,l] xsd;              // sample standard deviation of training obs

  int n_test;                   // number of test trials
  array[n_test] real x_test;          // locations of test trials
  array[n_test] int y_test;           // group labels for test trials
  array[n_test,m] int z_test_counts;  // responses for test trials
}

transformed data {
  real n_each;
  matrix[m,l] ss;
  ss = (n - 1) .* xsd;

  n_each = max(n);

}

parameters {
  // these are all shared across groups (same prior beliefs):
  real<lower=0> kappa_0;       // prior pseudocount for mean
  real<lower=0> nu_0;          // prior pseudocount for sd
  array[m] real m_0;                 // prior expected mean
  array[m] real<lower=0> S_0;        // prior expected standard deviation
  real<lower=0, upper=1> lapse_rate;
}

transformed parameters {
  // updated beliefs depend on input/group
  array[m,l] real m_n;                 // updated expected mean
  array[m,l] real<lower=0> kappa_n;    // updated mean pseudocount
  array[m,l] real<lower=0> S_n;        // updated expected sd
  array[m,l] real<lower=0> nu_n;       // updated sd pseudocount
  array[m,l] real<lower=0> t_scale;    // scale parameter of predictive t distribution
  array[n_test] simplex[m] p_test_conj;
  array[n_test] vector[m] log_p_test_conj;

  // update NIX2 parameters according to conjuate updating rules are taken from
  // Murphy (2007)
  for (cat in 1:m) {
    for (group in 1:l) {
      kappa_n[cat,group] = kappa_0 + n[cat,group];
      nu_n[cat,group] = nu_0 + n[cat,group];
      m_n[cat,group] = (m_0[cat] * kappa_0 + xbar[cat,group] * n[cat,group]) / kappa_n[cat,group];
      S_n[cat,group] = sqrt((nu_0*S_0[cat]^2 +
                                 ss[cat,group] +
                                 (n[cat,group]*kappa_0)/(kappa_n[cat,group]) *
                                   (m_0[cat] - xbar[cat,group])^2
                                 ) /
                                nu_n[cat,group]);
      t_scale[cat,group] = S_n[cat,group] * sqrt((kappa_n[cat,group] + 1) / kappa_n[cat,group]);
    }
  }

  // compute category probabilities for each of the test stimuli
  for (j in 1:n_test) {
    int group;
    group = y_test[j];
    // calculate un-normalized log prob for each category
    for (cat in 1:m) {
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
  vector[m] lapsing_probs;

  lapsing_probs = rep_vector(lapse_rate / m, m);

  // need to calculate category probabilities for each test trial
  kappa_0 ~ normal(0, n_each*4);
  nu_0 ~ normal(0, n_each*4);

  m_0 ~ normal(0, 100);
  S_0 ~ uniform(0, 100);

  for (i in 1:n_test) {
    z_test_counts[i] ~ multinomial(p_test_conj[i] * (1-lapse_rate) + lapsing_probs);
  }

}

generated quantities {
  // simplex[m] p_test_semi_conj[n_test];
  // simplex[m] p_test_sampled[n_test];
  // real ll;

  // ll = 0;
  // for (i in 1:n_test) {
  //   ll = ll + multinomial_log(z_test_counts[i], p_test_conj[i]);
  // }

}
