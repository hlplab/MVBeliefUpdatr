/*
  Fit multinomial response data using a belief-updating model to infer prior
  parameters.  A normal-inverse-Chisquare (NIX) prior is used, and it's assumed
  that the subject knows the true labels of all the input stimuli (e.g.,
  doesn't model any uncertainty in classification.

  Assumptions:

  (1) Beliefs about each category are assumed to have the form of multiple
      independent univariate NIX distributions, with unknown
      (inferred) m, S, kappa, and nu parameters describing uncertainty about the
      independent distributions of each cue for this category. During categorization,
      these independent beliefs about each cue are integrated (cue integration).
  (2) Beliefs for each category are independent of each other.
  (3) Subjects know the true labels of all the exposure stimuli.
  (4) Categorization has uniform decision-biases with unknown (inferred) lapse rate.

  Sufficient statistics and (if provided) prior mu_0 and Sigma_0 are assumed
  to be in the space defined by an affine transformation f(x) = SCALE(x + shift).
  Typically, that's the space obtained by whitening or standardizing the exposure
  data.

  Modified by
  Florian Jaeger
  April 2025
*/

data {
  int K;                                                   // number of features
  int L;                                                   // number of exposure groups (e.g. subjects)
  int M;                                                   // number of categories

  array[M,L] int<lower=0> N_exposure;                      // number of observations per category (M) and exposure group (L)
  array[M,L] vector[K] x_mean_exposure;                    // means for each category (M) and exposure group (L)
  array[M,L] cov_matrix[K] x_cov_exposure;                 // covariance matrix for each category (M) and group (L)

  int N_test;                                              // number of unique combinations of test locations & exposure groups
  array[N_test] vector[K] x_test;                          // locations (in cue space) of test trials
  array[N_test] int y_test;                                // exposure group indicator of test trials
  array[N_test,M] int z_test_counts;                       // responses for test trials (for each category M)

  simplex[M] phi;                                          // prior category probabilities

  int<lower=0, upper=1> lapse_rate_known;
  int<lower=0, upper=1> mu_0_known;
  int<lower=0, upper=1> sigma_0_known;
  /* Potential extension in the future: allow user-provided (weights_known = 1 & user_ideal_weights = 0) or
     inferred (weights_known = 0 & user_ideal_weights = 0) cue weights. Some lines below that are commented
     out contain drafts of this idea. Note, however, that those weights would be fixed (unless an updating
     mechanism is specified). So for now, we're assuming ideal weights, which do change with changes in the
     informativity of the cue.
  */
  // int<lower=0, upper=1> weights_known;
  // int<lower=0, upper=1> use_ideal_weights;                 // 0 = user-provided weights, 1 = ideal
  array[lapse_rate_known ? 1 : 0] real<lower=0, upper=1> lapse_rate_data;  // optional: user-provided lapse_rate
  array[mu_0_known ? M : 0] vector[mu_0_known ? K : 0] mu_0_data;          // optional: user-provided expected mu (prior category means) in space of affine transformation defined by INV_SCALE^-1 and shift
  array[sigma_0_known ? M : 0] vector[sigma_0_known ? K : 0] sigma_0_data; // optional: user-provided expected SD (prior category SD) in space of affine transformation defined by INV_SCALE^-1 and shift
  // simplex[use_ideal_weights || !weights_known ? 0 : K] cue_weights_data;   // optional: user-provided cue weights (used only if use_ideal_weights == 0 and weights are known)

  vector<lower=0>[mu_0_known ? 0 : K] tau_scale;           // separate taus for each of the K features to capture that features can be on separate scales

  /* The data above are assumed to have been transformed by a sensible affine transformation f(x) = A(x - center)
     (e.g., by scaling or whitening the data) in order to improve numerical stability. A_inv and center are used
     below in order to transform the parameters of interest back into the original space.
  */
  matrix[K, K] INV_SCALE;                                  // Inverse of transformation matrix (must be invertible) that was applied to data
  vector[K] shift;                                         // Center of affine transformation that was applied to data

  int<lower=0, upper=1> split_loglik_per_observation;
}

transformed data {
  array[M,L] vector[K] x_ss_exposure;                      // sum of *centered* squares matrix for each category (M) and group (L)
  for (category in 1:M) {
    for (group in 1:L) {
      for (cue in 1:K) {
        x_ss_exposure[category, group, cue] = x_cov_exposure[category, group, cue, cue] * (N_exposure[category, group] - 1);
      }
    }
  }

  real<lower=0> sigma_kappanu = max(max(to_array_1d(N_exposure)) * 4, 10); // scale for the prior of kappa/nu_0 (at least 10)
  vector[K] m_0_mu = rep_vector(0, K);                     // center of prior of m_0
}

parameters {
  /* The strength of the prior beliefs is currently assumed to be shared across groups, categories, and cues.
     Any of these assumptions could be revisited in the future. For example, different groups could differ in
     the strength of their prior beliefs because of differences in prior experience. Different categories could
     be associated with weaker or stronger prior beliefs because of differences in the amount of exposure
     (overall frequency of the category) or because differences in e.g., the cross-talker variability in the
     realization of the category. Similar reasoning can be applied to cues.
  */
  real<lower=K> kappa_0;                                        // prior pseudocount for category mu
  real<lower=K+1> nu_0;                                         // prior pseudocount for category Sigma

  // simplex[use_ideal_weights || weights_known ? 0 : K] cue_weights_param;         // cue weights parameter (used only if use_ideal_weights == 0)
  array[lapse_rate_known ? 0 : 1] real<lower=0, upper=1> lapse_rate_param;

  array[mu_0_known ? 0 : M] vector[K] m_0_param;                // prior expected means
  vector<lower=0>[mu_0_known ? 0 : K] m_0_tau;                  // prior SD of m_0

  array[sigma_0_known ? 0 : M] vector<lower=0>[K] s_0_param;    // square root of prior scale of Inverse Chisquare, s_0
}

transformed parameters {
  array[M] vector[K] m_0 = mu_0_known ? mu_0_data : m_0_param;                     // prior expected means
  array[M] vector<lower=0>[K] s_0;                                                 // square root of prior scale of Inverse Chisquare, s_0
  // simplex[use_ideal_weights ? 0 : K] cue_weights = weights_known ? cue_weights_data : cue_weights_param;
  real lapse_rate = lapse_rate_known ? lapse_rate_data[1] : lapse_rate_param[1];   // lapse rate

   /* Assuming unifom bias, so that lapsing_prob = probability of each category prior to
     taking into account stimulus is 1/M
  */
  vector[M] lapsing_probs = rep_vector(lapse_rate / M, M);

  // updated beliefs depend on input and group
  array[M,L] real<lower=K> kappa_n;          // updated mean pseudocount
  array[M,L] real<lower=K> nu_n;             // updated sd pseudocount
  array[M,L] vector[K] m_n;                  // updated m_0
  array[M,L] vector[K] s_n;                  // square root of updated scale of Inverse Chisquare, s_0
  array[M,L] vector[K] t_scale;              // scale matrix of predictive t distribution

  array[L] simplex[K] cue_weights;           // separate cue weights for each group (since the informativity of cues could differ between groups)

  array[N_test] simplex[M] p_test_conj;
  array[N_test] vector[M] log_p_test_conj;

  /* Update NIW parameters according to conjugate updating rules are taken from
     Murphy (2007, p. 136)
  */
  for (cat in 1:M) {
    // Get s_0 from expected Sigma given nu_0
    // E[sigma_c^2] = nu  / (nu - 2) * sigma_0^2 <==> sigma_0^2 = (E[sigma_c^2] * (nu_0 - 2)) / nu_0
    s_0[cat] = sigma_0_known ? sigma_0_data[cat] * sqrt((nu_0 - 2) / nu_0): s_0_param[cat];

    for (group in 1:L) {
      if (N_exposure[cat,group] > 0 ) {
        kappa_n[cat,group] = kappa_0 + N_exposure[cat,group];
        nu_n[cat,group] = nu_0 + N_exposure[cat,group];
        m_n[cat,group] =
          (kappa_0 * m_0[cat] + N_exposure[cat,group] * x_mean_exposure[cat,group]) /
          kappa_n[cat,group];
        s_n[cat,group] =
          sqrt(
             (nu_0 * s_0[cat]^2 +
              x_ss_exposure[cat,group] +
              (N_exposure[cat,group] * kappa_0) / (kappa_n[cat,group]) * (m_0[cat] - x_mean_exposure[cat,group])^2
          ) / nu_n[cat,group]
        );
      } else {
        kappa_n[cat,group] = kappa_0;
        nu_n[cat,group] = nu_0;
        m_n[cat,group] = m_0[cat];
        s_n[cat,group] = s_0[cat];
      }

      t_scale[cat,group] =
        s_n[cat,group] * sqrt((kappa_n[cat,group] + 1) / kappa_n[cat,group]);
    }
  }

  // Compute cue weights for each group
  for (group in 1:L) {
    vector[K] raw_weights = rep_vector(0, K);
    for (cat1 in 1:M) {
      for (cat2 in 1:M) {
        // NEEDS REVIEW: ARE IDEAL WEIGHTS REALLY A FUNCTION OF SIGMA_N (RATHER THAN, E.G., THE EXPECTED CATEGORY SD)?
        // Using element-wise operations ./ and .* so that this can all be done without looping over K
        raw_weights += phi[cat1] * phi[cat2] * ((m_n[cat1,group] - m_n[cat2,group])^2 ./ (s_n[cat1,group] .* s_n[cat2,group]));
      }
    }
    raw_weights ./= rep_vector(2.0, K);
    cue_weights[group] = softmax(raw_weights); // Enforces simplex
  }



  // compute posterior probabilities of all categories for each of the test stimuli
  for (j in 1:N_test) {
    int group;
    group = y_test[j];
    /* calculate un-normalized log posterior probability for each category. Under the assumption
       of uniform prior probabilities for each category, the log probabilities identical to the
       normalized log likelihoods. If we ever were to change this assumption, we'd have to add
       the log prior probabilities of categories here.
    */
    for (cat in 1:M) {
      // log likelihood of the test stimulus given the category
      log_p_test_conj[j,cat] = student_t_lpdf(x_test[j] |
                                              nu_n[cat,group],
                                              m_n[cat,group],
                                              t_scale[cat,group]);
      // NEEDS REVIEW: DOES THE CUE-WEIGHTING HAVE TO BE APPLIED BEFORE CALCULATING THE DENSITY?
      log_p_test_conj[j,cat] = sum(log_p_test_conj[j,cat] * cue_weights[group]);
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

      for (cat in 1:M) {
        m_0_param[cat] ~ normal(m_0_mu[cat], m_0_tau);
      }
  }

  // Specifying prior for components of S_0:
  if (!sigma_0_known) {
    for (cat in 1:M) {
      s_0_param[cat] ~ cauchy(0, tau_scale);
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

  // Inverse transform parameters back into the original space
  array[M] vector[K] m_0_original;
  array[M] vector[K] S_0_original;
  array[M,L] vector[K] m_n_original;
  array[M,L] vector[K] S_n_original;
  for (cat in 1:M) {
    m_0_original[cat] = INV_SCALE * m_0[cat] - shift;
    // NEEDS REVISITING: IS THIS REALLY THE CORRECT WAY TO BACKTRANSFORM A VECTOR OF (ESSENTIALLY) STANDARD DEVIATIONS?
    S_0_original[cat] = diagonal(INV_SCALE * diag_matrix(s_0[cat]^2) * INV_SCALE');               // store variance, not standard deviation (for the sake of parallelism to the NIW models)
    for (group in 1:L) {
      m_n_original[cat,group] = INV_SCALE * m_n[cat,group] - shift;
      // NEEDS REVISITING: IS THIS REALLY THE CORRECT WAY TO BACKTRANSFORM A VECTOR OF (ESSENTIALLY) STANDARD DEVIATIONS?
      S_n_original[cat,group] = diagonal(INV_SCALE * diag_matrix(s_n[cat,group]^2) * INV_SCALE'); // store variance, not standard deviation (for the sake of parallelism to the NIW models)
    }
  }
}
