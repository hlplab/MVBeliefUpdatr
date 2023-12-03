#' Get likelihood
#'
#' Get likelihood of observation(s) x given the exemplar model.
#'
#' @param x Observations. Can be a vector with k elements for a single observation or a matrix with k
#' columns and n rows, in which case each row of the matrix is taken to be one observation. If x is a
#' tibble with k columns or a list of vectors of length k, it is reduced into a matrix with k columns.
#' @param noise_treatment Determines whether perceptual noise is considered during categorization, and how.
#' Can be "no_noise", "sample", or "marginalize". If "no_noise", no noise will be applied to the input,
#' and no noise will be assumed during categorization. If "marginalize", average noise (i.e., no noise)
#' will be added to the stimulus, and `Sigma_noise` is added to Sigma when calculating the likelihood.
#' This simulates the expected consequences for perceptual noise on categorization *in the limit*, i.e,
#' if the input was categorized infinitely many times. If "sample", then noise is sampled and applied to
#' the input, and `Sigma_noise` is added to Sigma when calculating the likelihood. This simulates the
#' consequence of perceptual noise *on a particular observation*. If "sample" or "marginalize" are chosen,
#' `Sigma_noise` must be a covariance matrix of appropriate dimensions. (default: "no_noise" if Sigma_noise
#' is NULL, "marginalize" otherwise).
#' @param log Should the log-transformed density be returned (`TRUE`)? (default: `TRUE`)
#'
#' @seealso TBD
#' @keywords TBD
#' @importFrom rlang :=
#' @importFrom tidyr as_tibble
#' @export
get_likelihood_from_exemplars <- function(
  x,
  model,
  noise_treatment = if (is.exemplar_model(model)) { if (!is.null(first(model$Sigma_noise))) "marginalize" else "no_noise" } else "no_noise",
  log = T,
  category = "category",
  category.label = NULL
) {

}


#' Get categorization from an exemplar model
#'
#' Categorize a single observation based on an exemplar model. The decision rule can be specified to be either the
#' criterion choice rule, proportional matching (Luce's choice rule), or the sampling-based interpretation of
#' Luce's choice rule.
#'
#' @param x A vector of observations.
#' @param model An \code{\link[=is.exemplar_model]{exemplar_model}} object.
#' @param decision_rule Must be one of "criterion", "proportional", or "sampling".
#' @param noise_treatment Determines whether and how multivariate Gaussian noise is considered during categorization.
#' See \code{\link[=get_likelihood_from_exemplars]{get_likelihood_from_exemplars}}. (default: "sample" if decision_rule is "sample"; "marginalize" otherwise).
#' @param lapse_treatment Determines whether and how lapses will be treated. Can be "no_lapses", "sample" or "marginalize".
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
#' @export
get_categorization_from_exemplar_model <- function(
  x,
  model,
  decision_rule,
  noise_treatment = if (decision_rule == "sampling") "sample" else "marginalize",
  lapse_treatment = if (decision_rule == "sampling") "sample" else "marginalize",
  simplify = F
) {
  assert_that(is.exemplar_model(model))
  assert_that(decision_rule  %in% c("criterion", "proportional", "sampling"),
              msg = "Decision rule must be one of: 'criterion', 'proportional', or 'sampling'.")
  assert_that(any(lapse_treatment %in% c("no_lapses", "sample", "marginalize")),
              msg = "lapse_treatment must be one of 'no_lapses', 'sample' or 'marginalize'.")

  # In case a single x is handed as argument, make sure it's made a list so that the length check below
  # correctly treats it as length 1 (rather than the dimensionality of the one observation).
  if (!is.list(x)) x <- list(x)

  posterior_probabilities <-
    get_likelihood_from_exemplars(x = x, model = model, log = F, noise_treatment = noise_treatment) %>%
    group_by(category) %>%
    mutate(
      observationID = 1:length(x),
      x = x,
      lapse_rate = get_lapse_rate_from_model(model),
      lapse_bias = get_lapse_biases_from_model(model, categories = category),
      prior = get_priors_from_model(model, categories = category)) %>%
    group_by(observationID) %>%
    mutate(posterior_probability = (likelihood * prior) / sum(likelihood * prior))

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
    select(-c(likelihood, posterior_probability)) %>%
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
