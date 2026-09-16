# Update normalization -----------------------------------------------------
update_normalization_and_normalize_test_data <- function(
    mu_0 = first(prior_marginal_VOT_f0_stats$x_mean),
    kappa,
    data
) {
  # Get normalization parameters for each exposure data (exposure condition, list, block)
  exposure.normalization <-
    data %>%
    filter(Phase == "exposure") %>%
    group_by(ExposureGroup) %>%
    filter(ParticipantID == first(ParticipantID)) %>%
    summarise(
      x_N = length(x),
      x_mean = list(colMeans(reduce(x, rbind))),
      x_cov = list(cov(reduce(x, rbind)))) %>%
    # Apply normalization based on exposure to test
    mutate(
      mu_n = map2(x_N, x_mean, ~ 1 / (.env$kappa + .x) * (.env$kappa * .env$mu_0 + .x * .y))) %>%
    bind_rows(tibble(ExposureGroup = "no exposure", x_N = 0, x_mean = list(.env$mu_0), x_cov = NULL, mu_n = list(.env$mu_0)))

  # Apply normalization to each test data, and return those test data
  data %>%
    filter(Phase == "test") %>%
    left_join(
      exposure.normalization %>%
        select(ExposureGroup, mu_n),
      by = "ExposureGroup") %>%
    mutate(x = map2(x, mu_n, ~ .x - (.y - .env$mu_0)))
}

get_likelihood_from_updated_normalization <- function(
    par,
    prior = m_IO.VOT_f0,
    data = d_for_ASP.for_normalization
) {
  kappa <- exp(par[1])

  ll <-
    update_normalization_and_normalize_test_data(
      kappa = kappa,
      data = data) %>%
    group_by(ExposureGroup) %>%
    select(ExposureGroup, x, Response.Category) %>%
    nest(data = c(x, Response.Category)) %>%
    # Since the model never changes (only the cues do), there is only
    # one group of data for normalization
    mutate(model = list(prior)) %>%
    get_likelihood_from_grouped_data()

  history.optimization_normalization <<-
    bind_rows(
      history.optimization_normalization,
      tibble(kappa = kappa, log_likelihood = ll))

  return(ll)
}

# Update decision biases ---------------------------------------------------
#
# This function assumes that the error signal is relative to the expected response.
# This assumes that listeners can figure out---even on unlabelled trials---what the
# expected response is.
update_decision_bias_by_group <- function(
    prior,
    beta,
    data
) {
  cues <- get_cue_labels_from_model(prior)
  u <-
    data %>%
    group_map(
      .f = ~ update_model_decision_bias_incrementally(
        model = prior,
        beta = beta,
        exposure = .x,
        exposure.category = "Item.ExpectedResponse.Category",
        exposure.cues = cues,
        noise_treatment = "marginalize",
        lapse_treatment = "marginalize",
        keep.update_history = TRUE,
        keep.exposure_data = FALSE) %>%
        nest(posterior = everything()) %>%
        bind_cols(.y)) %>%
    reduce(bind_rows)
}

format_updated_bias_models_and_join_test_data <- function(
    models,
    prior,
    data.test
) {
  models %>%
    mutate(
      posterior =
        map(
          posterior,
          ~ filter(.x, observation.n %in% c(48, 96, 144)) %>%
            mutate(observation.n = case_when(
              observation.n == 48 ~ "_Up to test3",
              observation.n == 96 ~ "_Up to test5",
              observation.n == 144 ~ "_Up to test7")))) %>%
    unnest(posterior) %>%
    # Remap the different update steps onto the ExposureGroup
    mutate(ExposureGroup = paste0(gsub("^(.*)_.*$", "\\1", ExposureGroup), observation.n)) %>%
    select(-observation.n) %>%
    # If test data contains tests for "no exposure" condition, bind prior model back in
    # and call it "no exposure" (necessary for likelihood calculation)
    { if ("no exposure" %in% data.test$ExposureGroup)
      bind_rows(
        .,
        bind_cols(
          prior,
          tibble(ExposureGroup = "no exposure")) %>%
          crossing(ParticipantID = unique(as.character(data.test$ParticipantID)))) else . } %>%
    # Nest updated model and join test responses
    nest(model = -c(ExposureGroup, ParticipantID)) %>%
    left_join(data.test, by = join_by(ExposureGroup, ParticipantID))
}

get_likelihood_from_updated_bias <- function(
    par,
    prior = m_IO.VOT_f0,
    # These data sets are set as a default to speed up computation, since this allows us
    # to move any steps that would have to repeated on each optimization step but that do
    # *not* depend on beta out of the optimization.
    data.exposure = d_for_ASP.for_decision_changes.exposure,
    data.test = d_for_ASP.for_decision_changes.test
) {
  beta <- exp(par[1])

  ll <-
    suppressWarnings(
      update_decision_bias(
        prior = prior,
        beta = beta,
        data = data.exposure)) %>%
    format_updated_bias_models_and_join_test_data(prior = prior, data.test = data.test) %>%
    get_likelihood_from_grouped_data()

  message("For beta = ", beta, " found log-likelihood of ", ll)
  history.optimization_bias <<-
    bind_rows(
      history.optimization_bias,
      tibble(
        beta = beta,
        log_likelihood = ll))

  return(ll)
}

# Make trivariate normal distribution ---------------------------------------------------
# Function takes as arguments the mean and covariance matrix of a trivariate normal distribution over (in this order)
# VOT, f0_Mel, and vowel_duration, as well as values for f0_Mel and vowel_duration. The function then returns a
# *function* that represents the conditional normal distribution over VOT at that combination of f0_Mel and
# vowel_duration values. This function takes VOT as input and returns a density value.
#
# This function was drafted with the help of Google's Bard: https://bard.google.com/chat/32a2428ef7e806bb
conditional_univariate_normal_from_trivariate_normal <- function(mu, Sigma, f0_Mel, vowel_duration) {
  # Calculate the conditional mean and conditional standard deviation of x given y = f0_Mel and z = vowel_duration
  mu_x <- mu[1] + Sigma[1, 2] / Sigma[2, 2] * (f0_Mel - mu[2]) + Sigma[1, 3] / Sigma[3, 3] * (vowel_duration - mu[3])
  sigma_x <- sqrt(Sigma[1, 1] - Sigma[1, 2] / Sigma[2, 2] * Sigma[2, 1] - Sigma[1, 3] / Sigma[3, 3] * Sigma[3, 1])

  # Define the function for the conditional univariate normal distribution over x given f0_Mel and vowel_duration
  f <- function(x) {
    1 / (sqrt(2 * pi) * sigma_x) * exp(-0.5 * ((x - mu_x) / sigma_x)^2)
  }

  return(f)
}

# This functions takes as arguments the means and covariance matrices of two trivariate normal distributions,
# as well as values for f0_Mel and vowel_duration. The function then returns a *function* that represents the
# posterior of the category corresponding to the second trivariate normal distribution.  This function takes
# VOT as input and returns a posterior probability.
conditional_univariate_posterior_t <- function(mu_d, Sigma_d, mu_t, Sigma_t, f0_Mel, vowel_duration) {
  density_t <- conditional_univariate_normal_from_trivariate_normal(mu_t, Sigma_t, f0_Mel, vowel_duration)
  density_d <- conditional_univariate_normal_from_trivariate_normal(mu_d, Sigma_d, f0_Mel, vowel_duration)

  f <- function(x) {
    density_t(x) / (density_t(x) + density_d(x))
  }
}
