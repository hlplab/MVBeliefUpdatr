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
  noise_treatment = infer_default_noise_treatment(model$Sigma_noise),
  log = T,
  category = "category",
  category.label = NULL
) {

}

# Deprecated after S7-migration

#' Legacy wrapper for categorize.
#'
#' @description Deprecated. Use \code{\link{categorize}} instead.
#' @rdname get_categorization_from_model
#' @export
#' @description Deprecated. Use categorize() instead.
#' @keywords internal
get_categorization_from_exemplar_model <- function(
  x,
  model,
  decision_rule,
  noise_treatment = if (decision_rule == "sampling") "sample" else infer_default_noise_treatment(model$Sigma_noise),
  lapse_treatment = if (decision_rule == "sampling") "sample" else "marginalize",
  simplify = F
) {
  lifecycle::deprecate_warn(
    when = "0.0.3",
    what = "get_categorization_from_exemplar_model()",
    with = "categorize()"
  )
  .assert_that(is.exemplar_model(model))
  .assert_that(decision_rule  %in% c("criterion", "proportional", "sampling"),
              msg = "Decision rule must be one of: 'criterion', 'proportional', or 'sampling'.")
  .assert_that(any(lapse_treatment %in% c("no_lapses", "sample", "marginalize")),
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
      lapse_rate = get_lapse_rate(model),
      lapse_bias = get_lapse_bias(model, categories = category),
      prior = get_category_prior(model, categories = category)) %>%
    group_by(observationID) %>%
    mutate(posterior_probability = (likelihood * prior) / sum(likelihood * prior))

  # How should lapses be treated?
  if (lapse_treatment == "sample") {
    posterior_probabilities %<>%
      mutate(
        posterior_probability = ifelse(
          rep(
            rbinom(1, 1, lapse_rate),
            length(get_category_labels(model))),
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
            length(get_category_labels(model))),
          posterior_probability + runif(
            length(get_category_labels(model)),
            min = 0,
            max = 0),
          posterior_probability),
        # select most probable category
        response = ifelse(posterior_probability == max(posterior_probability), 1, 0))
  } else if (decision_rule == "sampling") {
    posterior_probabilities %<>%
      mutate(response = .rmultinom(1, 1, posterior_probability) %>% as.vector())
  } else if (decision_rule == "proportional") {
    posterior_probabilities %<>%
      mutate(response = posterior_probability)
  } else warning("Unsupported decision rule. This should be impossible to happen. Do not trust the results.")

  posterior_probabilities %<>%
    ungroup() %>%
    select(-c(likelihood, posterior_probability)) %>%
    select(observationID, x, category, response)

  if (simplify) {
    .assert_that(decision_rule  %in% c("criterion", "sampling"),
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
