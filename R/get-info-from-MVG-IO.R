#' Get likelihood
#'
#' Get likelihood of observation(s) x given the MVG parameters mu and Sigma. This is the density of
#' a multivariate normal distribution over k dimensions.
#'
#' @param x Observations. Can be a vector with k elements for a single observation or a matrix with k
#' columns and n rows, in which case each row of the matrix is taken to be one observation. If x is a
#' tibble with k columns or a list of vectors of length k, it is reduced into a matrix with k columns.
#' @param mu The category mean mu. Should be a k x 1 or 1 x k
#' matrix, or vector of length k.
#' @param Sigma The category covariance matrix Sigma. Should be a square k x k matrix.
#' @param Sigma_noise Optionally, a covariance matrix describing the perceptual noise to be applied while
#' calculating the posterior predictive. (default: `NULL`)
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
#' @rdname get_MVG_likelihood
#' @export
get_MVG_likelihood <- function(
    x, mu, Sigma, Sigma_noise = NULL,
    noise_treatment = infer_default_noise_treatment(Sigma_noise),
    log = T
) {
  # mvtnorm::dmvt expects means to be vectors, and x to be either a vector or a matrix.
  # in the latter case, each *row* of the matrix is an input.
  .assert_that(is.vector(x) | is.matrix(x) | is_tibble(x) | is.list(x))
  .assert_that(is.vector(mu) | is.matrix(mu) | .is_non_NA_scalar_double(mu))
  .assert_that(.is_sigma(Sigma))

  # do not reorder these conditionals (go from more to less specific)
  if (is.matrix(mu)) mu <- as.vector(mu)
  x %<>% format_input_for_likelihood_calculation(dim = length(mu))
  .assert_that(dim(x)[2] == length(mu),
              msg = "Input x and m are not of compatible dimensions.")

  .assert_that(.is_non_NA_scalar_logical(log))
  .assert_that(any(noise_treatment %in% c("no_noise", "marginalize", "sample")),
              msg = "noise_treatment must be one of 'no_noise', 'marginalize', or 'sample'.")
  if (noise_treatment != "no_noise") {
    .assert_that(.is_sigma(Sigma_noise),
                msg = 'If noise_treatment is not "no_noise", Sigma_noise must be a covariance matrix of appropriate dimensions, matching those of the category covariance matrices Sigma.')
    .assert_that(all(dim(Sigma) == dim(Sigma_noise)),
                msg = 'If noise_treatment is not "no_noise", Sigma_noise must be a covariance matrix of appropriate dimensions, matching those of the category covariance matrices Sigma.')
  }

  D <- get_D(Sigma)
  if (D == 1) {
    .assert_that(.is_non_NA_scalar_double(mu), msg = "Sigma and mu are not of compatible dimensions.")
  } else {
    .assert_that(dim(Sigma)[2] == D,
                msg = "Sigma is not a square matrix, and thus not a covariance matrix")
    .assert_that(length(mu) == dim(x)[2],
                msg = paste("mu and input are not of compatible dimensions. mu is of length", length(mu), "but input has", dim(x)[2], "columns."))
    .assert_that(length(mu) == D,
                msg = "Sigma and mu are not of compatible dimensions.")
  }

  if (noise_treatment == "sample") {
    .assert_that(
      is_weakly_greater_than(nrow(x), 1),
      msg = "For noise sampling, x must be of length 1 or longer.")

    x <- x + .rmvnorm(n = nrow(x), mean = rep(0, ncol(x)), sigma = Sigma_noise)
  }

  if (noise_treatment %in% c("sample", "marginalize")) {
    Sigma <- Sigma + Sigma_noise
  }

  .dmvnorm(x, mean = mu, sigma = Sigma, log = log) %>% as.numeric()
}


#' @rdname get_MVG_likelihood
#' @importFrom rlang :=
#' @importFrom tidyr pivot_wider as_tibble
#' @export
get_likelihood_from_MVG <- function(
  x,
  model,
  noise_treatment = infer_default_noise_treatment(model$Sigma_noise),
  log = T,
  category = "category",
  category.label = NULL,
  wide = FALSE
) {
  .assert_that(is.MVG(model))
  .assert_optional_character(category.label)
  .assert_that(any(noise_treatment == "no_noise", is.MVG_ideal_observer(model)),
              msg = 'No noise matrix Sigma_noise found. If noise_treatment is not "no_noise", then model must be an MVG_ideal_observer.')

  if (is.null(category.label)) {
    model %<>%
      droplevels()

    category.label <-
      model %>%
      dplyr::pull(!! sym(category)) %>%
      unique()
  }

  likelihood <- foreach(c = category.label) %do% {
    m <-
      model %>%
      filter(!! sym(category) == c)

    get_MVG_likelihood(
      x = x,
      mu = m$mu[[1]],
      Sigma = m$Sigma[[1]],
      log = log,
      noise_treatment = noise_treatment,
      Sigma_noise = if (noise_treatment == "no_noise") NULL else m$Sigma_noise[[1]]) %>%
      as_tibble(.name_repair = "unique") %>%
      rename_with(~ if (log) { "log_likelihood" } else { "likelihood" }) %>%
      mutate(!! sym(category) := c)
  }
  likelihood %<>% reduce(rbind)

  if (wide)
    likelihood %<>%
    pivot_wider(
      values_from = if (log) "log_likelihood" else "likelihood",
      names_from = !! sym(category),
      names_prefix = if (log) "log_likelihood." else "likelihood.") %>%
    unnest()

  return(likelihood)
}

# Deprecated after S7-migration

#' Legacy wrapper for posterior.
#'
#' @description Deprecated. Use \code{\link{posterior}} instead.
#' @rdname get_posterior_from_model
#' @export
#' @description Deprecated. Use posterior() instead.
#' @keywords internal
get_posterior_from_MVG_ideal_observer <- function(
    x,
    model,
    noise_treatment = if (decision_rule == "sampling") "sample" else infer_default_noise_treatment(model$Sigma_noise),
    lapse_treatment = if (decision_rule == "sampling") "sample" else "marginalize"
) {
  lifecycle::deprecate_warn(
    when = "0.0.3",
    what = "get_posterior_from_MVG_ideal_observer()",
    with = "posterior()",
    always = TRUE
  )

  # TO DO: check dimensionality of x with regard to belief.
  assert_MVG_ideal_observer(model)
  .assert_that(any(lapse_treatment %in% c("no_lapses", "sample", "marginalize")),
              msg = "lapse_treatment must be one of 'no_lapses', 'sample' or 'marginalize'.")

  # When the input isn't a list, that's ambiguous between the input being a single input or a set of
  # 1D inputs. Use the model's cue dimensionality to disambiguate between the two cases.
  if (!is.list(x)) {
    x <- if (get_cue_dimensionality_from_model(model) == 1) as.list(x) else list(x)
  } else if (
    length(x) > 0L &&
    all(vapply(x, function(value) length(value) == 1L && is.atomic(value), logical(1)))
  ) {
    x <- list(unlist(x, use.names = FALSE))
  }

  x_input <- x
  n.distinct_categories <- length(get_category_labels(model))
  if (!is.list(x_input)) {
    x_input <- if (get_cue_dimensionality_from_model(model) == 1) as.list(x_input) else list(x_input)
  }

  posterior_probabilities <-
    get_likelihood_from_MVG(x = x, model = model, log = F, noise_treatment = noise_treatment) %>%
    tibble::as_tibble() %>%
    mutate(
      observationID = rep(seq_along(x_input), times = n.distinct_categories),
      x = rep(x_input, times = n.distinct_categories)
    )

  lapse_rate <- get_lapse_rate(model)
  posterior_probabilities$lapse_rate <- lapse_rate
  posterior_probabilities$lapse_bias <- get_lapse_bias(model, categories = posterior_probabilities$category)
  posterior_probabilities$prior <- get_category_prior(model, categories = posterior_probabilities$category)

  posterior_probabilities <-
    posterior_probabilities %>%
    group_by(observationID) %>%
    mutate(posterior_probability = (.data$likelihood * .data$prior) / sum(.data$likelihood * .data$prior))

  # How should lapses be treated?
  if (lapse_treatment == "sample") {
    posterior_probabilities %<>%
      mutate(
        posterior_probability = ifelse(
          rep(
            rbinom(1, 1, lapse_rate),
            length(get_category_labels(model))),
          .data$lapse_bias,                 # substitute lapse probabilities for posterior
          .data$posterior_probability))     # ... or not
  } else if (lapse_treatment == "marginalize") {
    posterior_probabilities %<>%
      mutate(posterior_probability = lapse_rate * .data$lapse_bias + (1 - lapse_rate) * .data$posterior_probability)
  }

  posterior_probabilities %<>%
    ungroup() %>%
    select(-c(likelihood)) %>%
    select(observationID, x, category, posterior_probability) %>%
    arrange(.data$observationID)

  # Warn if any posteriors don't sum up to 1.
  posterior.check <-
    posterior_probabilities %>%
    group_by(x, observationID) %>%
    summarise(posterior_probability = sum(.data$posterior_probability))

  posterior.check %<>%
    arrange(posterior_probability) %>%
    filter(is.na(posterior_probability) | is.nan(posterior_probability) | !isTRUE(all.equal(posterior_probability, 1.0, tolerance = 1e-8)))

  if (nrow(posterior.check) > 0L) {
    s <- paste(
      nrow(posterior.check),
      "input(s) have an ill-defined posterior under the model. This can happen when inputs are far away from all category means.\n")
    posterior.check %<>%
      mutate(
        string = purrr::pmap_chr(
          .l = list(posterior_probability, x, observationID),
          .f = function(posterior_probability, x, observationID) {
            paste0("Sum of posterior is ", posterior_probability, " for observation ID = ", observationID, "; input = ", paste(x, collapse = ","))
          }
        )
      )
    s %<>% paste0(., paste(posterior.check$string, collapse = ".\n"))
    warning(s)
  }

  return(posterior_probabilities)
}

#' Legacy wrapper for categorize.
#'
#' @description Deprecated. Use \code{\link{categorize}} instead.
#' @rdname get_categorization_from_model
#' @export
#' @description Deprecated. Use categorize() instead.
#' @keywords internal
get_categorization_from_MVG_ideal_observer <- function(
  x,
  model,
  decision_rule = "sampling",
  noise_treatment = if (decision_rule == "sampling") "sample" else infer_default_noise_treatment(model$Sigma_noise),
  lapse_treatment = if (decision_rule == "sampling") "sample" else "marginalize",
  simplify = F
) {
  lifecycle::deprecate_warn(
    when = "0.0.3",
    what = "get_categorization_from_MVG_ideal_observer()",
    with = "categorize()",
    always = TRUE
  )

  posterior_probabilities <-
    get_posterior_from_MVG_ideal_observer(x = x, model = model, noise_treatment = noise_treatment, lapse_treatment = lapse_treatment)

  # Apply decision rule
  if (decision_rule == "criterion") {
    posterior_probabilities %<>%
      group_by(observationID, x) %>%
      mutate(
        # tie breaker in case of uniform probabilities
        posterior_probability = ifelse(
          rep(
            sum(.data$posterior_probability == max(.data$posterior_probability)) > 1,
            length(get_category_labels(model))),
          posterior_probability + runif(
            length(get_category_labels(model)),
            min = 0,
            max = 1),
          .data$posterior_probability),
        # select most probable category
        response = ifelse(.data$posterior_probability == max(.data$posterior_probability), 1, 0))
  } else if (decision_rule == "sampling") {
    posterior_probabilities %<>%
      group_by(observationID, x) %>%
      mutate(response = .rmultinom(1, 1, .data$posterior_probability) %>% as.vector())
  } else if (decision_rule == "proportional") {
    posterior_probabilities %<>%
      mutate(response = .data$posterior_probability)
  } else warning("Unsupported decision rule. This should be impossible to happen. Do not trust the results.")

  posterior_probabilities %<>%
    ungroup() %>%
    select(-c(posterior_probability)) %>%
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
