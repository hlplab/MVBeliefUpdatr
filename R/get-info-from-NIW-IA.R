#' @rdname get_categorization_function
#' @export
get_categorization_function_from_NIW_ideal_adaptor = function(x, ...) {
  get_categorization_function(
    ms = x$m,
    Ss = x$S,
    kappas = x$kappa,
    nus = x$nu,
    lapse_rate = x$lapse_rate,
    lapse_bias = x$lapse_bias,
    ...
  )
}

#' Get categorization from an NIW ideal adaptor
#'
#' Categorize a single observation based on an NIW ideal adaptor The decision rule can be specified to be either the
#' criterion choice rule, proportional matching (Luce's choice rule), or the sampling-based interpretation of
#' Luce's choice rule.
#'
#' @param x A vector of observations.
#' @param model An \code{\link[=is.NIW_ideal_adaptor]{NIW_ideal_adaptor}} object.
#' @param decision_rule Must be one of "criterion", "proportional", or "sampling".
#' @param noise_treatment Determines how multivariate Gaussian noise is added to the input. Can be "no_noise", "sample"
#' or "marginalize". If "sample", observations are adjusted by samples drawn from the noise distribution before applying
#' categorization.If "marginalize" then each observation is transformed into the marginal distribution
#' that results from convolving the input with noise. This latter option might be helpful, for example, if one is
#' interested in estimating the consequences of noise across individuals. (default: "sample" if decision_rule is
#' "sample"; "marginalize" otherwise).
#' @param lapse_treatment Determines how lapses will be treated. Can be "no_lapses", "sample" or "marginalize".
#' If "sample", whether a trial is lapsing or not will be sampled for each observations. If a trial is sampled to be
#' a lapsing trial the lapse biases are used as the posterior for that trial. If "marginalize", the posterior probability
#' will be adjusted based on the lapse formula lapse_rate * lapse_bias + (1 - lapse_rate) * posterior probability from
#' perceptual model. (default: "sample" if decision_rule is "sample"; "marginalize" otherwise).
#' @param simplify Should the output be simplified, and just the label of the selected category be returned? This
#' option is only available for the criterion and sampling decision rules. (default: `FALSE`)
#'
#' @return Either a tibble of observations with posterior probabilities for each category (in long format), or a
#' character vector indicating the chosen category in the same order as the observations in x (if simplify = `TRUE`).
#'
#' @seealso TBD
#' @keywords TBD
#' @examples
#' TBD
#' @export
#'
get_categorization_from_NIW_ideal_adaptor = function(
  x,
  model,
  decision_rule,
  noise_treatment = if (decision_rule == "sampling") "sample" else "marginalize",
  lapse_treatment = if (decision_rule == "sampling") "sample" else "marginalize",
  simplify = F
) {
  # TO DO: check dimensionality of x with regard to model.
  assert_NIW_ideal_adaptor(model)
  assert_that(decision_rule  %in% c("criterion", "proportional", "sampling"),
              msg = "Decision rule must be one of: 'criterion', 'proportional', or 'sampling'.")
  assert_that(any(noise_treatment %in% c("no_noise", "sample", "marginalize")),
              msg = "noise_treatment must be one of 'no_noise', 'sample' or 'marginalize'.")
  assert_that(any(lapse_treatment %in% c("no_noise", "sample", "marginalize")),
              msg = "lapse_treatment must be one of 'no_noise', 'sample' or 'marginalize'.")

  # In case a single x is handed as argument, make sure it's made a list so that the length check below
  # correctly treats it as length 1 (rather than the dimensionality of the one observation).
  if (!is.list(x)) x <- list(x)

  # How should noise be treated?
  if (noise_treatment == "sample") {
    assert_that(
      is_weakly_greater_than(length(x), 1),
      msg = "For noise sampling, x must be of length 1 or longer.")

    x <- map(x, ~ rmvnorm(n = 1, mean = .x, sigma = model$Sigma_noise[[1]]))
  } else if (noise_treatment == "marginalize") {
      model %<>%
        mutate(Sigma = map2(Sigma, Sigma_noise, ~ .x + .y))
  }

  posterior_probabilities <-
    get_posterior_predictive_from_NIW_belief(x = x, model = model, log = F) %>%
    group_by(category) %>%
    mutate(
      observationID = 1:length(x),
      x = x,
      lapse_rate = get_lapse_rate_from_model(model),
      lapse_bias = get_lapse_biases_from_model(model, categories = category),
      prior = get_priors_from_model(model, categories = category)) %>%
    group_by(observationID) %>%
    mutate(posterior_probability = (posterior_predictive * prior) / sum(posterior_predictive * prior))

  # How should lapses be treated?
  if (lapse_treatment == "sample") {
    posterior_probabilities %<>%
      mutate(
        posterior_probability = ifelse(
          rep(
            rbinom(1, 1, lapse_rate),
            get_nlevels_of_category_labels_from_model(model)),
          lapse_bias,                 # substitute lapse probabilities for posterior
          posterior_probability))     # ... or not
  } else if (lapse_treatment == "marginalize") {
    posterior_probabilities %<>%
      mutate(posterior_probability =  lapse_rate * lapse_bias + (1 - lapse_rate) * posterior_probability)
  }

  # Apply decision rule
  if (decision_rule == "criterion") {
    posterior_probabilities %<>%
      mutate(
        # tie breaker in case of uniform probabilities
        posterior_probability = ifelse(
          rep(
            sum(posterior_probability == max(posterior_probability)) > 1,
            get_nlevels_of_category_labels_from_model(model)),
          posterior_probability + runif(
            get_nlevels_of_category_labels_from_model(model),
            min = 0,
            max = 0),
          posterior_probability),
        # select most probable category
        response = ifelse(posterior_probability == max(posterior_probability), 1, 0))
  } else if (decision_rule == "sampling") {
    posterior_probabilities %<>%
      mutate(response = rmultinom(1, 1, posterior_probability) %>% as.vector())
  } else if (decision_rule == "proportional") {
    posterior_probabilities %<>%
      mutate(response = posterior_probability)
  } else warning("Unsupported decision rule. This should be impossible to happen. Do not trust the results.")

  posterior_probabilities %<>%
    ungroup() %>%
    select(-c(posterior_predictive, posterior_probability)) %>%
    select(observationID, x, category, response)

  if (simplify) {
    assert_that(decision_rule  %in% c("criterion", "sampling"),
                msg = "For simplify = T, decision rule must be either criterion or sampling.")
    return(posterior_probabilities %>%
             filter(response == 1) %>%
             select(observationID, category) %>%
             arrange(observationID) %>%
             rename(response = category) %>%
             ungroup() %>%
             pull(response))
  } else return(posterior_probabilities)
}
