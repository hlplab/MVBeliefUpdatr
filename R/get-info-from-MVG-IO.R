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
  noise_treatment = "no_noise",
  log = T,
  category = "category",
  category.label = NULL,
  wide = FALSE
) {
  .assert_optional_character(category.label)

  d <- length(get_cue_labels(model))
  x_mat <- .legacy_observation_matrix(x, d)
  lik <- likelihood(model, x_mat, categories = category.label)
  if (log) lik <- log(lik)

  cats <- colnames(lik)
  value_name <- if (log) "log_likelihood" else "likelihood"

  if (wide) {
    out <- tibble::as_tibble(as.data.frame(lik), .name_repair = "minimal")
    names(out) <- paste0(value_name, ".", cats)
    return(out)
  }

  out <- tibble(
    !! sym(value_name) := as.vector(lik),
    !! sym(category) := rep(cats, each = nrow(lik))
  )
  out
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
    decision_rule = "sampling",
    noise_treatment = if (decision_rule == "sampling") "sample" else "no_noise",
    lapse_treatment = if (decision_rule == "sampling") "sample" else "marginalize"
) {
  lifecycle::deprecate_warn(
    when = "0.0.3",
    what = "get_posterior_from_MVG_ideal_observer()",
    with = "posterior()",
    always = TRUE
  )

  .legacy_long_posterior(model, x, noise_treatment, lapse_treatment)
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
  noise_treatment = if (decision_rule == "sampling") "sample" else "no_noise",
  lapse_treatment = if (decision_rule == "sampling") "sample" else "marginalize",
  simplify = F
) {
  lifecycle::deprecate_warn(
    when = "0.0.3",
    what = "get_categorization_from_MVG_ideal_observer()",
    with = "categorize()",
    always = TRUE
  )

  d.response <-
    .legacy_long_posterior(model, x, noise_treatment, lapse_treatment) %>%
    .legacy_apply_decision_rule(decision_rule)

  if (simplify) .legacy_simplify_categorization(d.response, decision_rule) else d.response
}
