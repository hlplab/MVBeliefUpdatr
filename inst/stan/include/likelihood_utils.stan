/*
  Likelihood and Mixture Utility Functions for MVBeliefUpdatr
  
  Covers:
  - Mixture of perceptual category decision rule and lapse response
  - Split observation-level log-likelihood computation
*/

// Categorization response probability under uniform/biased lapse mixture
vector lapse_mixture(vector p_choice, real lapse_rate, vector lapse_bias) {
  return p_choice * (1.0 - lapse_rate) + lapse_rate * lapse_bias;
}

/* Compute observation-level split log-likelihood for LOO/WAIC diagnostics. This is
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

     Optionally, users can thus also request log-likelihoods to be calculated for each individual response.
     This requires reformatting the information z_test_counts. The downside of this revised approach is that
     it makes the object very large to store, and slow to work with.
  */
vector compute_split_loglik(
  array[,] int z_counts,
  array[] vector p_choice,
  real lapse_rate,
  vector lapse_bias,
  int split_loglik,
  int total_obs
) {
  if (split_loglik == 0) {
    return rep_vector(0.0, 0);
  }
  vector[total_obs] loglik_split;
  int idx = 1;
  int n_test = size(p_choice);
  int M = num_elements(lapse_bias);
  for (n in 1:n_test) {
    vector[M] p_mix = lapse_mixture(p_choice[n], lapse_rate, lapse_bias);
    for (m in 1:M) {
      real log_p = log(p_mix[m]);
      for (i in 1:z_counts[n, m]) {
        loglik_split[idx] = log_p;
        idx += 1;
      }
    }
  }
  return loglik_split;
}
