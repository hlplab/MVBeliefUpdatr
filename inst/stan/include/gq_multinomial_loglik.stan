  /* Compute and store pointwise log-likelihood for LOOIC / WAIC computation */
  real log_lik_sum = 0;
  vector[n_test] log_lik;

  for (n in 1:n_test) {
    log_lik[n] = multinomial_lpmf(z_test_counts[n] | lapse_mixture(p_test_conj[n], lapse_rate, lapse_bias));
    log_lik_sum += log_lik[n];
  }

  vector[split_loglik_per_observation ? sum(to_array_1d(z_test_counts)) : 0] log_lik_split =
    compute_split_loglik(
      z_test_counts,
      p_test_conj,
      lapse_rate,
      lapse_bias,
      split_loglik_per_observation,
      split_loglik_per_observation ? sum(to_array_1d(z_test_counts)) : 0
    );
