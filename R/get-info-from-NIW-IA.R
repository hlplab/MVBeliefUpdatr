#' Get NIW categorization function
#'
#' Returns a categorization function for the first category, based on a set of parameters for the Normal-Inverse-Wishart (NIW)
#' distribution. ms, Ss, kappas, nus, and priors are assumed to be of the same length and sorted the same way, so that the first
#' element of ms is corresponding to the same category as the first element of Ss, kappas, nus, and priors, etc.
#'
#' @param ms Means of the multivariate normal distributions over category means.
#' @param Ss Scatter matrices of the inverse Wishart distribution over category covariance matrices.
#' @param kappas Strength of the beliefs into the distribution over category means.
#' @param nus Strength of the beliefs into the distribution over category covariance matrices.
#' @param priors Vector of categories' prior probabilities. (default: uniform prior over categories)
#' @param lapse_rate A lapse rate for the categorization responses.
#' @param lapse_biases A lapse bias for the categorization responses. (default: uniform bias over categories)
#' @param Sigma_noise A noise matrix. (default: a 0-matrix)
#' @param noise_treatment How should the noise specified in \code{Sigma_noise} be considered in the categorization function?
#' For details, see \code{\link{get_NIW_posterior_predictive}}.
#' @param logit Should the function that is returned return log-odds (TRUE) or probabilities (FALSE)? (default: TRUE)
#'
#' @return A function that takes as input cue values and returns posterior probabilities of the first category,
#' based on the posterior predictive of the cues given the (IBBU-derived parameters for the) categories' m, S,
#' kappa, nu, and prior, as well as the lapse rate.
#'
#' @seealso TBD
#' @keywords TBD
#'
#' @rdname get_NIW_categorization_function
#' @export
get_NIW_categorization_function = function(
    ms, Ss, kappas, nus,
    priors = rep(1 / length(ms), length(ms)),
    lapse_rate = NULL,
    lapse_biases = rep(1 / length(ms), length(ms)),
    Sigma_noise = NULL,
    noise_treatment = if (any(is.null(Sigma_noise), all(is.null(Sigma_noise)), all(map_lgl(Sigma_noise, is.null)))) "no_noise" else "marginalize",
    logit = FALSE
) {
  tolerance = 1e-5
  assert_that(are_equal(length(ms), length(Ss)),
              are_equal(length(ms), length(priors)),
              are_equal(length(ms), length(kappas)),
              are_equal(length(ms), length(nus)),
              msg = "The number of ms, Ss, kappas, nus, and priors must be identical.")
  n.cat = length(ms)

  assert_that(all(between(priors, 0, 1), between(sum(priors), 1 - tolerance, 1 + tolerance)),
              msg = "priors must sum to 1.")
  if (any(is.null(lapse_rate),
          all(is.null(lapse_rate)),
          all(map_lgl(lapse_rate, is.null)))) {
    lapse_rate = rep(0, n.cat)
  } else {
    assert_that(all(between(lapse_rate, 0, 1)))
  }
  if (any(is.null(lapse_biases),
          all(is.null(lapse_biases)),
          all(map_lgl(lapse_biases, is.null)))) {
    lapse_biases <- 1 / n.cat
  } else {
    assert_that(all(between(lapse_biases, 0, 1), between(sum(lapse_biases), 1 - tolerance, 1 + tolerance)),
                msg = "lapse biases must sum to 1.")
  }
  if (any(is.null(Sigma_noise),
          all(is.null(Sigma_noise)),
          all(map_lgl(Sigma_noise, is.null)))) {
    Sigma_noise <-
      matrix(
        0,
        nrow = if (is.null(dim(Ss[[1]]))) 1 else max(dim(Ss[[1]])),
        ncol = if (is.null(dim(Ss[[1]]))) 1 else max(dim(Ss[[1]])))
  } else {
    assert_that(is.Sigma(Sigma_noise))
  }

  # Get dimensions of multivariate category
  D = get_D(ms)
  assert_that(
    nus[[1]] >= D,
    msg = "Nu must be at least K (number of dimensions of the multivariate Gaussian category).")

  f <- function(x, target_category = 1) {
    if (!is.list(x)) x <- list(x)
    log_p <- matrix(nrow = length(x), ncol = n.cat) # this seems to assume that x is a list
    for (cat in 1:n.cat) {
      log_p[, cat] <-
        get_NIW_posterior_predictive(
          x, # can this handle lists?
          ms[[cat]], Ss[[cat]], kappas[[cat]], nus[[cat]],
          Sigma_noise = Sigma_noise, noise_treatment = noise_treatment,
          log = T)
    }

    p_target <-
      (1 - lapse_rate[target_category]) * exp(log_p[,target_category] + log(priors[target_category]) - log(rowSums(exp(log_p) * priors))) +
      lapse_rate[target_category] * lapse_biases[target_category]

    if (logit)
      return(qlogis(p_target))
    else
      return(p_target)
  }

  return(f)
}


#' @rdname get_NIW_categorization_function
#' @export
get_categorization_function_from_NIW_ideal_adaptor = function(model, ...) {
  get_NIW_categorization_function(
    ms = model$m,
    Ss = model$S,
    kappas = model$kappa,
    nus = model$nu,
    priors = model$prior,
    lapse_rate = model$lapse_rate,
    lapse_biases = model$lapse_bias,
    Sigma_noise = model$Sigma_noise,
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
#' @param noise_treatment Determines whether and how multivariate Gaussian noise is considered during categorization.
#' See \code{\link[=get_NIW_posterior_predictive]{get_NIW_posterior_predictive}}. (default: "sample" if decision_rule is
#' "sample"; "marginalize" otherwise).
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
#'
get_categorization_from_NIW_ideal_adaptor = function(
  x,
  model,
  decision_rule,
  noise_treatment = if (decision_rule == "sampling") "sample" else "marginalize",
  lapse_treatment = if (decision_rule == "sampling") "sample" else "marginalize",
  simplify = F,
  verbose = F
) {
  # TO DO: check dimensionality of x with regard to model.
  assert_NIW_ideal_adaptor(model, verbose = verbose)
  assert_that(decision_rule  %in% c("criterion", "proportional", "sampling"),
              msg = "Decision rule must be one of: 'criterion', 'proportional', or 'sampling'.")
  assert_that(any(lapse_treatment %in% c("no_lapses", "sample", "marginalize")),
              msg = "lapse_treatment must be one of 'no_lapses', 'sample' or 'marginalize'.")

  # In case a single x is handed as argument, make sure it's made a list so that the length check below
  # correctly treats it as length 1 (rather than the dimensionality of the one observation).
  if (!is.list(x)) x <- list(x)

  n.distinct_categories <- length(get_category_labels_from_model(model))
  posterior_probabilities <-
    get_posterior_predictive_from_NIW_belief(x = x, model = model, log = F, noise_treatment = noise_treatment) %>%
    mutate(
      observationID = rep(1:length(.env$x), .env$n.distinct_categories),
      x = rep(.env$x, .env$n.distinct_categories),
      lapse_rate = get_lapse_rate_from_model(.env$model),
      lapse_bias = get_lapse_biases_from_model(.env$model, categories = .data$category),
      prior = get_priors_from_model(.env$model, categories = .data$category)) %>%
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
