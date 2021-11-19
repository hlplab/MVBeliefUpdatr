#' Get likelihood
#'
#' Get likelihood of observations x given the MVG parameters mu and Sigma. This is the density of
#' a multivariate normal distribution.
#'
#' @param x Observations. Can be a vector with k elements for a single observation or a matrix with k
#' columns and n rows, in which case each row of the matrix is taken to be one observation. If x is a
#' tibble with k columns or a list of vectors of length k, it is reduced into a matrix with k columns.
#' @param mu The category mean mu. Should be a
#' matrix or vector of length k.
#' @param Sigma The category covariance matrix Sigma. Should be a square k x k matrix.
#' @param log Should the log-transformed density be returned (`TRUE`)? (default: `TRUE`)
#' @param noise_treatment Determines whether and how multivariate Gaussian noise is added to the input. Can be "no_noise", "sample"
#' or "marginalize". If "no_noise", no noise will be applied. If "sample" or "marginalize", `Sigma_noise` must be a covariance
#' matrix of appropriate dimensions. If "sample", observations are adjusted by samples drawn from the noise distribution before applying
#' categorization. If "marginalize" then each observation is transformed into the marginal distribution
#' that results from convolving the input with noise. This latter option might be helpful, for example, if one is
#' interested in estimating the consequences of noise across individuals. (default: "no_noise").
#' @param Sigma_noise Optionally, a covariance matrix describing the perceptual noise to be applied while
#' calculating the posterior predictive. (default: `NULL`)
#'
#' @seealso TBD
#' @keywords TBD
#' @examples
#' TBD
#' @rdname get_MVG_likelihood
#' @export
get_MVG_likelihood = function(x, mu, Sigma, log = T, noise_treatment = "no_noise", Sigma_noise = NULL) {
  # mvtnorm::dmvt expects means to be vectors, and x to be either a vector or a matrix.
  # in the latter case, each *row* of the matrix is an input.
  assert_that(is.vector(x) | is.matrix(x) | is_tibble(x) | is.list(x))
  assert_that(is.vector(mu) | is.matrix(mu) | is_scalar_double(mu))
  assert_that(is.positive.definite(Sigma) | is_scalar_double(Sigma))
  # do not reorder these conditionals (go from more to less specific)
  if (is_tibble(x)) x %<>% as.matrix() else
    if (is.list(x)) x %<>% reduce(rbind) %<>% as.matrix(nrow = 1) else
      if (is.vector(x)) x %<>% matrix(nrow = 1)
  if (is.matrix(mu)) mu = as.vector(mu)

  assert_that(is.flag(log))
  assert_that(any(noise_treatment %in% c("no_noise", "sample", "marginalize")),
              msg = "noise_treatment must be one of 'no_noise', 'sample' or 'marginalize'.")
  if (noise_treatment != "no_noise") {
    assert_that(is.positive.definite(Sigma_noise)  | is_scalar_double(Sigma_noise))
    assert_that(all(dim(Sigma) == dim(Sigma_noise)),
                msg = 'If noise_treatment is not "no_noise", Sigma_noise must be a covariance matrix of appropriate dimensions, matching those of the category covariance matrices Sigma.')
  }

  D = dim(Sigma)[1]

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

  # How should noise be treated?
  if (noise_treatment == "sample") {
    assert_that(
      is_weakly_greater_than(length(x), 1),
      msg = "For noise sampling, x must be of length 1 or longer.")

    x <- map(x, ~ rmvnorm(n = 1, mean = .x, sigma = Sigma_noise))
  } else if (noise_treatment == "marginalize") {
    Sigma = Sigma + Sigma_noise
  }

  dmvnorm(x,
       mean = mu,
       sigma = Sigma,
       log = log)
}


#' @rdname get_MVG_likelihood
#' @export
get_likelihood_from_MVG = function(
  x,
  model,
  log = T,
  noise_treatment = if (is.MVG_ideal_adaptor(model)) "marginalize" else "no_noise",
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
    category.label = model %>%
      pull(!! sym(category)) %>%
      unique()
  }

  likelihood <- foreach(c = category.label) %do% {
    m <-
      model %>%
      filter(!! sym(category) == c)

    get_likelihood(
      x = x,
      mu = m$mu[[1]],
      Sigma = m$Sigma[[1]],
      log = log,
      noise_treatment = noise_treatment,
      Sigma_noise = if (noise_treatment == "no_noise") NULL else m$Sigma_noise[[1]]) %>%
      as_tibble() %>%
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


#' Get categorization from an ideal observer
#'
#' Categorize a single observation based on an MVG ideal observer The decision rule can be specified to be either the
#' criterion choice rule, proportional matching (Luce's choice rule), or the sampling-based interpretation of
#' Luce's choice rule.
#'
#' @param x A vector of observations.
#' @param model An \code{\link[=is.MVG_ideal_observer]{MVG_ideal_observer}} object.
#' @param decision_rule Must be one of "criterion", "proportional", or "sampling".
#' @param noise_treatment Determines whether and how multivariate Gaussian noise is added to the input.
#' See \code{\link[=get_MVG_likelihood]{get_MVG_likelihood}}. (default: "sample" if decision_rule is
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
#' @examples
#' TBD
#' @export
get_categorization_from_MVG_ideal_observer = function(
  x,
  model,
  decision_rule,
  noise_treatment = if (decision_rule == "sampling") "sample" else "marginalize",
  lapse_treatment = if (decision_rule == "sampling") "sample" else "marginalize",
  simplify = F
) {
  # TO DO: check dimensionality of x with regard to belief.
  assert_MVG_ideal_observer(model)
  assert_that(decision_rule  %in% c("criterion", "proportional", "sampling"),
              msg = "Decision rule must be one of: 'criterion', 'proportional', or 'sampling'.")
  assert_that(any(lapse_treatment %in% c("no_noise", "sample", "marginalize")),
              msg = "lapse_treatment must be one of 'no_noise', 'sample' or 'marginalize'.")

  # In case a single x is handed as argument, make sure it's made a list so that the length check below
  # correctly treats it as length 1 (rather than the dimensionality of the one observation).
  if (!is.list(x)) x <- list(x)

  posterior_probabilities <-
    get_likelihood_from_MVG(x = x, model = model, log = F, noise_treatment = noise_treatment) %>%
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
