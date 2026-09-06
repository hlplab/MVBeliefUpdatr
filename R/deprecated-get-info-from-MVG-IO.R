#' @include S7-core-classes.R
#' @importFrom lifecycle deprecate_warn
#' @importFrom rlang := sym
#' @importFrom tidyr pivot_wider as_tibble
#' @importFrom mvtnorm rmvnorm
NULL

# deprecated ------------------------------------------------------------------

#' Deprecated: get_MVG_likelihood
#'
#' @description `r lifecycle::badge("deprecated")`
#' `get_MVG_likelihood()` was deprecated in MVBeliefUpdatr 0.1.0 and will be
#' removed in 0.2.0. Please use [likelihood()] instead.
#'
#' @param x Observations.
#' @param mu Mean vector.
#' @param Sigma Covariance matrix.
#' @param Sigma_noise Optional noise covariance matrix.
#' @param noise_treatment Noise treatment (`"no_noise"`, `"sample"`, or `"marginalize"`).
#' @param log Logical; whether log likelihood is returned.
#' @return Numeric vector of likelihoods.
#' @seealso [likelihood()]
#' @keywords internal
#' @rdname get_MVG_likelihood
#' @export
get_MVG_likelihood <- function(
  x,
  mu,
  Sigma,
  Sigma_noise = NULL,
  noise_treatment = .infer_noise_treatment(Sigma_noise),
  log = TRUE
) {
  lifecycle::deprecate_warn(
    "0.1.0",
    "get_MVG_likelihood()",
    with = "likelihood()"
  )

  .assert_that(is.vector(mu) || is.matrix(mu) || .is_scalar_double(mu))
  .assert_that(.is_sigma(Sigma))
  if (is.matrix(mu)) mu <- as.vector(mu)

  d <- length(mu)
  x <- .as_observation_matrix(x, d = d, arg_name = "x")
  .assert_that(.is_scalar_logical(log))
  .assert_that(any(noise_treatment %in% c("no_noise", "marginalize", "sample")),
    msg = "noise_treatment must be one of 'no_noise', 'marginalize', or 'sample'."
  )

  if (noise_treatment != "no_noise") {
    .assert_that(.is_sigma(Sigma_noise))
    .assert_that(all(dim(Sigma) == dim(Sigma_noise)))
  }

  if (noise_treatment == "sample") {
    .assert_that(nrow(x) >= 1, msg = "For noise sampling, x must be of length 1 or longer.")
    x <- x + mvtnorm::rmvnorm(n = nrow(x), mean = rep(0, ncol(x)), sigma = Sigma_noise)
  }

  if (noise_treatment %in% c("sample", "marginalize")) {
    Sigma <- Sigma + Sigma_noise
  }

  .dmvnorm(x, mean = mu, sigma = Sigma, log = log) %>% as.numeric()
}

#' Deprecated: get_likelihood_from_MVG
#'
#' @description `r lifecycle::badge("deprecated")`
#' `get_likelihood_from_MVG()` was deprecated in MVBeliefUpdatr 0.1.0 and will be
#' removed in 0.2.0. Please use [likelihood()] instead.
#'
#' @param x Observations.
#' @param model An MVG model.
#' @param noise_treatment Noise treatment.
#' @param log Logical; whether log likelihood is returned.
#' @param category Category column name.
#' @param category.label Category labels.
#' @param wide Logical; whether wide format is returned.
#' @return Likelihood data frame.
#' @seealso [likelihood()]
#' @keywords internal
#' @rdname get_MVG_likelihood
#' @export
get_likelihood_from_MVG <- function(
  x,
  model,
  noise_treatment = "no_noise",
  log = TRUE,
  category = "category",
  category.label = NULL,
  wide = FALSE
) {
  lifecycle::deprecate_warn(
    "0.1.0",
    "get_likelihood_from_MVG()",
    with = "likelihood()"
  )

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

  tibble::tibble(
    !!rlang::sym(value_name) := as.vector(lik),
    !!rlang::sym(category) := rep(cats, each = nrow(lik))
  )
}

#' Deprecated: get_posterior_from_MVG_ideal_observer
#'
#' @description `r lifecycle::badge("deprecated")`
#' `get_posterior_from_MVG_ideal_observer()` was deprecated in MVBeliefUpdatr 0.1.0 and will be
#' removed in 0.2.0. Please use [posterior()] instead.
#'
#' @param x Observations.
#' @param model Model object.
#' @param decision_rule Decision rule.
#' @param noise_treatment Noise treatment.
#' @param lapse_treatment Lapse treatment.
#' @return Posterior data frame.
#' @seealso [posterior()]
#' @keywords internal
#' @rdname get_posterior_from_model
#' @export
get_posterior_from_MVG_ideal_observer <- function(
  x,
  model,
  decision_rule = "sampling",
  noise_treatment = if (decision_rule == "sampling") "sample" else "no_noise",
  lapse_treatment = if (decision_rule == "sampling") "sample" else "marginalize"
) {
  lifecycle::deprecate_warn(
    "0.1.0",
    "get_posterior_from_MVG_ideal_observer()",
    with = "posterior()"
  )

  .legacy_long_posterior(model, x, noise_treatment, lapse_treatment)
}

#' Deprecated: get_categorization_from_MVG_ideal_observer
#'
#' @description `r lifecycle::badge("deprecated")`
#' `get_categorization_from_MVG_ideal_observer()` was deprecated in MVBeliefUpdatr 0.1.0 and will be
#' removed in 0.2.0. Please use [categorize()] instead.
#'
#' @param x Observations.
#' @param model Model object.
#' @param decision_rule Decision rule.
#' @param noise_treatment Noise treatment.
#' @param lapse_treatment Lapse treatment.
#' @param simplify Logical; whether to simplify to category vector.
#' @return Categorization data frame or vector.
#' @seealso [categorize()]
#' @keywords internal
#' @rdname get_categorization_from_model
#' @export
get_categorization_from_MVG_ideal_observer <- function(
  x,
  model,
  decision_rule = "sampling",
  noise_treatment = if (decision_rule == "sampling") "sample" else "no_noise",
  lapse_treatment = if (decision_rule == "sampling") "sample" else "marginalize",
  simplify = FALSE
) {
  lifecycle::deprecate_warn(
    "0.1.0",
    "get_categorization_from_MVG_ideal_observer()",
    with = "categorize()"
  )

  d.response <-
    .legacy_long_posterior(model, x, noise_treatment, lapse_treatment) %>%
    .legacy_apply_decision_rule(decision_rule)

  if (simplify) .legacy_simplify_categorization(d.response, decision_rule) else d.response
}

