#' Example MVG ideal observers.
#'
#' @export
example_MVG_ideal_observer <- function(example = 1) {
  if (example == 1) {
    message("An example MVG ideal observer for two categories in a 2D cue continuum that differ in means and correlatation, but not standard deviations. Lapse rate is .05 with uniform prior and lapse bias. No perceptual noise.")
    tibble(
      category = c("A", "B"),
      mu = list(c("cue1" = -2, "cue2" = -2), c("cue1" = 2, "cue2" = 2)),
      Sigma = list(
        matrix(c(3, 2.4, 2.4, 3), nrow = 2, dimnames = list(c("cue1", "cue2"), c("cue1", "cue2"))),
               matrix(c(3, -2.4, -2.4, 3), nrow = 2, dimnames = list(c("cue1", "cue2"), c("cue1", "cue2")))),
      prior = c(.5, .5),
      lapse_rate = .05,
      lapse_bias = c(.5, .5),
      Sigma_noise = list(
        matrix(rep(0, 4), nrow = 2, dimnames = list(c("cue1", "cue2"), c("cue1", "cue2"))),
        matrix(rep(0, 4), nrow = 2, dimnames = list(c("cue1", "cue2"), c("cue1", "cue2"))))) %>%
      mutate(across(category, factor))
  } else if (example == 2) {
    message("An example MVG ideal observer for two categories in a 1D cue continuum. Lapse rate is .05 with uniform prior and lapse bias. No perceptual noise.")
    tibble(
      category = c("/d/", "/t/"),
      mu = list(c("VOT" = 5), c("VOT" = 50)),
      Sigma = list(matrix(c(80), nrow = 1, dimnames = list(c("VOT"), c("VOT"))),
                   matrix(c(340), nrow = 1, dimnames = list(c("VOT"), c("VOT")))),
      prior = c(.5, .5),
      lapse_rate = .05,
      lapse_bias = c(.5, .5),
      Sigma_noise = list(
        matrix(c(0), nrow = 1, dimnames = list(c("VOT"), c("VOT"))),
        matrix(c(0), nrow = 1, dimnames = list(c("VOT"), c("VOT"))))) %>%
      mutate(category = factor(category))
  }
}


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
#' @examples
#' TBD
#' @rdname get_MVG_likelihood
#' @export
get_MVG_likelihood <- function(
    x, mu, Sigma, Sigma_noise = NULL,
    noise_treatment = if (is.null(Sigma_noise)) "no_noise" else "marginalize",
    log = T
) {
  # mvtnorm::dmvt expects means to be vectors, and x to be either a vector or a matrix.
  # in the latter case, each *row* of the matrix is an input.
  assert_that(is.vector(x) | is.matrix(x) | is_tibble(x) | is.list(x))
  assert_that(is.vector(mu) | is.matrix(mu) | is_scalar_double(mu))
  assert_that(is.Sigma(Sigma))

  # do not reorder these conditionals (go from more to less specific)
  if (is.matrix(mu)) mu = as.vector(mu)
  x %<>% format_input_for_likelihood_calculation()
  assert_that(dim(x)[2] == length(mu),
              msg = "Input x and m are not of compatible dimensions.")

  assert_that(is.flag(log))
  assert_that(any(noise_treatment %in% c("no_noise", "marginalize", "sample")),
              msg = "noise_treatment must be one of 'no_noise', 'marginalize', or 'sample'.")
  if (noise_treatment != "no_noise") {
    assert_that(is.Sigma(Sigma_noise),
                msg = 'If noise_treatment is not "no_noise", Sigma_noise must be a covariance matrix of appropriate dimensions, matching those of the category covariance matrices Sigma.')
    assert_that(all(dim(Sigma) == dim(Sigma_noise)),
                msg = 'If noise_treatment is not "no_noise", Sigma_noise must be a covariance matrix of appropriate dimensions, matching those of the category covariance matrices Sigma.')
  }

  D = get_D(Sigma)
  if (D == 1) {
    assert_that(is_scalar_double(mu), msg = "Sigma and mu are not of compatible dimensions.")
  } else {
    assert_that(dim(Sigma)[2] == D,
                msg = "Sigma is not a square matrix, and thus not a covariance matrix")
    assert_that(length(mu) == dim(x)[2],
                msg = paste("mu and input are not of compatible dimensions. mu is of length", length(mu), "but input has", dim(x)[2], "columns."))
    assert_that(length(mu) == D,
                msg = "Sigma and mu are not of compatible dimensions.")
  }

  if (noise_treatment == "sample") {
    assert_that(
      is_weakly_greater_than(nrow(x), 1),
      msg = "For noise sampling, x must be of length 1 or longer.")

    x <- x + rmvnorm(n = nrow(x), mean = rep(0, ncol(x)), sigma = Sigma_noise)
  }

  if (noise_treatment %in% c("sample", "marginalize")) {
    Sigma <- Sigma + Sigma_noise
  }

  dmvnorm(x, mean = mu, sigma = Sigma, log = log) %>% as.numeric()
}


#' @rdname get_MVG_likelihood
#' @importFrom rlang :=
#' @importFrom tidyr pivot_wider as_tibble
#' @export
get_likelihood_from_MVG <- function(
  x,
  model,
  noise_treatment = if (is.MVG_ideal_observer(model)) { if (!is.null(first(model$Sigma_noise))) "marginalize" else "no_noise" } else "no_noise",
  log = T,
  category = "category",
  category.label = NULL,
  wide = FALSE
) {
  assert_that(is.MVG(model))
  assert_that(any(is.null(category.label) | is.character(category.label)))
  assert_that(any(noise_treatment == "no_noise", is.MVG_ideal_observer(model)),
              msg = 'No noise matrix Sigma_noise found. If noise_treatment is not "no_noise", then model must be an MVG_ideal_observer.')

  if (is.null(category.label)) {
    model %<>%
      droplevels()

    category.label <-
      model %>%
      pull(!! sym(category)) %>%
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


#' @rdname get_posterior_from_model
#' @export
get_posterior_from_MVG_ideal_observer <- function(
    x,
    model,
    noise_treatment = if (decision_rule == "sampling") "sample" else "marginalize",
    lapse_treatment = if (decision_rule == "sampling") "sample" else "marginalize"
) {
  # TO DO: check dimensionality of x with regard to belief.
  assert_MVG_ideal_observer(model)
  assert_that(any(lapse_treatment %in% c("no_lapses", "sample", "marginalize")),
              msg = "lapse_treatment must be one of 'no_lapses', 'sample' or 'marginalize'.")

  # When the input isn't a list, that's ambiguous between the input being a single input or a set of
  # 1D inputs. Use the model's cue dimensionality to disambiguate between the two cases.
  if (!is.list(x)) {
    x <- if (get_cue_dimensionality_from_model(model) == 1) as.list(x) else list(x)
  }

  n.distinct_categories <- length(get_category_labels_from_model(model))
  posterior_probabilities <-
    get_likelihood_from_MVG(x = x, model = model, log = F, noise_treatment = noise_treatment) %>%
    mutate(
      observationID = rep(1:length(.env$x), .env$n.distinct_categories),
      x = rep(.env$x, .env$n.distinct_categories),
      lapse_rate = get_lapse_rate_from_model(.env$model),
      lapse_bias = get_lapse_biases_from_model(.env$model, categories = .data$category),
      prior = get_priors_from_model(.env$model, categories = .data$category)) %>%
    group_by(observationID) %>%
    mutate(posterior_probability = (.data$likelihood * .data$prior) / sum(.data$likelihood * .data$prior))

  # How should lapses be treated?
  if (lapse_treatment == "sample") {
    posterior_probabilities %<>%
      mutate(
        posterior_probability = ifelse(
          rep(
            rbinom(1, 1, .data$lapse_rate),
            get_nlevels_of_category_labels_from_model(model)),
          .data$lapse_bias,                 # substitute lapse probabilities for posterior
          .data$posterior_probability))     # ... or not
  } else if (lapse_treatment == "marginalize") {
    posterior_probabilities %<>%
      mutate(posterior_probability = .data$lapse_rate * .data$lapse_bias + (1 - .data$lapse_rate) * .data$posterior_probability)
  }

  posterior_probabilities %<>%
    ungroup() %>%
    select(-c(likelihood)) %>%
    select(observationID, x, category, posterior_probability)

  # Warn if any posteriors don't sum up  to 1.
  posterior.check <-
    posterior_probabilities %>%
    group_by(x, observationID) %>%
    summarise(posterior_probability = sum(.data$posterior_probability))
  if (any(
    is.na(posterior.check$posterior_probability),
    is.nan(posterior.check$posterior_probability),
    !near(posterior.check$posterior_probability, 1.0))) {
    posterior.check %<>%
      arrange(posterior_probability) %>%
      filter(is.na(posterior_probability) | is.nan(posterior_probability) | !near(posterior_probability, 1.0))
    s <- paste(
      nrow(posterior.check),
      "input(s) have an ill-defined posterior under the model. This can happen when inputs are far away from all category means.\n")
    s %<>% paste0(
      .,
      posterior.check %>%
        arrange(posterior_probability) %>%
        filter(is.na(posterior_probability) | is.nan(posterior_probability) | !near(posterior_probability, 1.0)) %>%
        mutate(string = pmap(
          .l = list(posterior_probability, x, observationID),
          .f = ~ paste0("Sum of posterior is ", ..1, " for observation ID = ", ..3, "; input = ", paste(..2, collapse = ","))) %>%
            unlist()) %>%
        pull(string) %>% paste(., collapse = ".\n"))
    warning(s)
  }

  return(posterior_probabilities)
}


#' @rdname get_categorization_from_model
#' @export
get_categorization_from_MVG_ideal_observer <- function(
  x,
  model,
  decision_rule,
  noise_treatment = if (decision_rule == "sampling") "sample" else "marginalize",
  lapse_treatment = if (decision_rule == "sampling") "sample" else "marginalize",
  simplify = F
) {
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
            get_nlevels_of_category_labels_from_model(model)),
          posterior_probability + runif(
            get_nlevels_of_category_labels_from_model(model),
            min = 0,
            max = 1),
          .data$posterior_probability),
        # select most probable category
        response = ifelse(.data$posterior_probability == max(.data$posterior_probability), 1, 0))
  } else if (decision_rule == "sampling") {
    posterior_probabilities %<>%
      group_by(observationID, x) %>%
      mutate(response = rmultinom(1, 1, .data$posterior_probability) %>% as.vector())
  } else if (decision_rule == "proportional") {
    posterior_probabilities %<>%
      mutate(response = .data$posterior_probability)
  } else warning("Unsupported decision rule. This should be impossible to happen. Do not trust the results.")

  posterior_probabilities %<>%
    ungroup() %>%
    select(-c(posterior_probability)) %>%
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
