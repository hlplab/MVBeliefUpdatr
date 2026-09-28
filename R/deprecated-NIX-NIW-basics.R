#' @include S7-class.R
#' @importFrom lifecycle deprecate_warn
#' @importFrom mvtnorm dmvt rmvnorm
NULL

# =============================================================================
# Deprecated NIX / NIW Functions
# =============================================================================

#' Deprecated: get_D
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#' `get_D()` is deprecated. Use `length(get_cue_labels(x))` instead.
#'
#' @param x A model, template, representation, or object from which cue labels can be extracted.
#' @return Number of cues (dimensionality) as an integer.
#' @keywords internal
#' @export
get_D <- function(x) {
  lifecycle::deprecate_warn(
    "0.1.0",
    "get_D()",
    details = "Use length(get_cue_labels(x)) instead."
  )
  length(get_cue_labels(x))
}

#' Deprecated: get_NIW_posterior_predictive
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#' Get posterior predictive density of observations x given the Normal-Inverse-Wishart (NIW)
#' parameters m, S, kappa, and nu. This is the density of a multivariate Student-T distribution
#' \insertCite{@see @murphy2012 p. 134}{MVBeliefUpdatr}.
#'
#' @param x Observation(s). Can be a vector with k elements, a matrix with k columns, a data frame,
#'   or a list of numeric vectors.
#' @param m Mean vector of length k.
#' @param S Scatter matrix of dimension k x k.
#' @param kappa Strength of belief on the mean (pseudocounts).
#' @param nu Strength of belief on the covariance matrix (pseudocounts).
#' @param Sigma_noise Optional perceptual noise covariance matrix (default: `NULL`).
#' @param noise_treatment Noise treatment: `"no_noise"`, `"sample"`, or `"marginalize"`.
#' @param log Logical; if `TRUE`, return log-transformed density (default: `TRUE`).
#'
#' @references \insertRef{murphy2012}{MVBeliefUpdatr}
#' @rdname get_NIW_posterior_predictive
#' @keywords internal
#' @export
get_NIW_posterior_predictive <- function(
  x,
  m,
  S,
  kappa,
  nu,
  Sigma_noise = NULL,
  noise_treatment = .infer_noise_treatment(Sigma_noise),
  log = TRUE
) {
  lifecycle::deprecate_warn(
    "0.1.0",
    "get_NIW_posterior_predictive()"
  )
  .assert_that(is.vector(m) || is.matrix(m) || .is_scalar_double(m))
  .assert_that(is.matrix(S) || .is_scalar_numeric(S))
  if (is.matrix(m)) m <- as.vector(m)

  d <- length(m)
  x <- .as_observation_matrix(x, d = d, arg_name = "x")

  .assert_that(all(.is_scalar_numeric(kappa), .is_scalar_numeric(nu)))
  .assert_that(.is_scalar_logical(log))
  .assert_that(any(noise_treatment %in% c("no_noise", "sample", "marginalize")),
    msg = "noise_treatment must be one of 'no_noise', 'sample' or 'marginalize'."
  )

  if (noise_treatment != "no_noise") {
    .assert_that(.is_sigma(Sigma_noise))
    .assert_that(all(dim(S) == dim(Sigma_noise)),
      msg = "Unless noise_treatment is 'no_noise', Sigma_noise must be a covariance matrix of appropriate dimensions."
    )
  }

  D <- if (is.matrix(S)) ncol(S) else length(m)
  .assert_that(nu >= D,
    msg = "nu must be at least as large as the number of dimensions of the multivariate Normal."
  )

  if (noise_treatment == "sample") {
    .assert_that(nrow(x) >= 1, msg = "For noise sampling, x must be of length 1 or longer.")
    x <- x + mvtnorm::rmvnorm(n = nrow(x), mean = rep(0, ncol(x)), sigma = Sigma_noise)
  }

  if (noise_treatment %in% c("sample", "marginalize")) {
    S <- get_S_from_expected_Sigma(get_expected_Sigma_from_S(S, nu) + Sigma_noise, nu)
  }

  scale_mat <- S * ((kappa + 1) / (kappa * (nu - D + 1)))
  df_val <- nu - D + 1

  .dmvt(x, delta = m, sigma = scale_mat, df = df_val, log = log)
}

#' Deprecated: get_NIX_posterior_predictive
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#' Get posterior predictive density of 1D observations x given the Normal-Inverse-Chi-Squared (NIX)
#' parameters m, sigma2 (or S), kappa, and nu. In 1D, this evaluates a Student-T density.
#'
#' @param x Observation(s) (numeric vector, matrix with 1 column, or 1D data frame).
#' @param m Prior mean location.
#' @param sigma2 Prior expected variance parameter (or variance scale \eqn{\sigma_0^2}).
#' @param kappa Strength of belief on the mean (pseudocounts).
#' @param nu Degrees of freedom / strength of belief on variance (pseudocounts).
#' @param Sigma_noise Optional 1x1 noise covariance / variance (default: `NULL`).
#' @param noise_treatment Noise treatment: `"no_noise"`, `"sample"`, or `"marginalize"`.
#' @param log Logical; if `TRUE`, return log-transformed density (default: `TRUE`).
#'
#' @rdname get_NIX_posterior_predictive
#' @keywords internal
#' @export
get_NIX_posterior_predictive <- function(
  x,
  m,
  sigma2,
  kappa,
  nu,
  Sigma_noise = NULL,
  noise_treatment = .infer_noise_treatment(Sigma_noise),
  log = TRUE
) {
  lifecycle::deprecate_warn(
    "0.1.0",
    "get_NIX_posterior_predictive()"
  )
  .assert_that(.is_scalar_numeric(m))
  .assert_that(.is_scalar_numeric(sigma2))
  .assert_that(.is_scalar_numeric(kappa))
  .assert_that(.is_scalar_numeric(nu))
  .assert_that(.is_scalar_logical(log))

  x <- .as_observation_matrix(x, d = 1, arg_name = "x")

  if (identical(noise_treatment, "sample") && !is.null(Sigma_noise)) {
    x <- x + mvtnorm::rmvnorm(n = nrow(x), mean = 0, sigma = as.matrix(Sigma_noise))
  }

  noise_variance <- if (!is.null(Sigma_noise) && (identical(noise_treatment, "sample") || identical(noise_treatment, "marginalize"))) {
    as.numeric(Sigma_noise)[1]
  } else {
    0
  }

  scale_eff <- sqrt((as.numeric(sigma2) * (as.numeric(kappa) + 1) / as.numeric(kappa)) + noise_variance)
  z <- (x[, 1] - as.numeric(m)) / scale_eff
  nu0 <- as.numeric(nu)

  if (isTRUE(log)) {
    stats::dt(z, df = nu0, log = TRUE) - log(scale_eff)
  } else {
    stats::dt(z, df = nu0) / scale_eff
  }
}
