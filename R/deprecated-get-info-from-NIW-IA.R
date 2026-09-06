#' @include S7-core-classes.R
#' @importFrom lifecycle deprecate_warn
#' @importFrom purrr map_lgl
NULL

# deprecated ------------------------------------------------------------------

#' Deprecated: get_NIW_categorization_function
#'
#' @description `r lifecycle::badge("deprecated")`
#' `get_NIW_categorization_function()` was deprecated in MVBeliefUpdatr 0.1.0 and will be
#' removed in 0.2.0. Please use [get_category_posterior_function()] instead.
#'
#' @param ms Means of the multivariate normal distributions over category means.
#' @param Ss Scatter matrices of the inverse Wishart distribution over category covariance matrices.
#' @param kappas Strength of the beliefs into the distribution over category means.
#' @param nus Strength of the beliefs into the distribution over category covariance matrices.
#' @param priors Vector of categories' prior probabilities.
#' @param lapse_rate A lapse rate for the categorization responses.
#' @param lapse_biases A lapse bias for the categorization responses.
#' @param Sigma_noise A noise matrix.
#' @param noise_treatment Noise treatment.
#' @param lapse_treatment Lapse treatment.
#' @return A categorization function.
#' @seealso [get_category_posterior_function()]
#' @keywords internal
#' @rdname get_NIW_categorization_function
#' @export
get_NIW_categorization_function <- function(
  ms, Ss, kappas, nus,
  priors = rep(1 / length(ms), length(ms)),
  lapse_rate = 0,
  lapse_biases = rep(1 / length(ms), length(ms)),
  Sigma_noise = NULL,
  noise_treatment = .infer_noise_treatment(Sigma_noise),
  lapse_treatment = if (lapse_rate > 0) "marginalize" else "no_lapses"
) {
  lifecycle::deprecate_warn(
    "0.1.0",
    "get_NIW_categorization_function()",
    with = "get_category_posterior_function()"
  )

  tolerance <- MVBU_PROB_TOL
  .assert_that(.is_equal(length(ms), length(Ss)),
    .is_equal(length(ms), length(priors)),
    .is_equal(length(ms), length(kappas)),
    .is_equal(length(ms), length(nus)),
    msg = "The number of ms, Ss, kappas, nus, and priors must be identical."
  )
  n.cat <- length(ms)

  .assert_that(all(between(priors, 0, 1), between(sum(priors), 1 - tolerance, 1 + tolerance)),
    msg = "priors must sum to 1."
  )
  .assert_that(.is_scalar_double(lapse_rate),
    msg = "lapse_rate must be a scalar."
  )
  .assert_that(between(lapse_rate, 0, 1))
  if (any(is.null(lapse_biases),
    all(is.null(lapse_biases)),
    all(purrr::map_lgl(lapse_biases, is.null)))) {
    lapse_biases <- 1 / n.cat
  } else {
    .assert_that(all(between(lapse_biases, 0, 1), between(sum(lapse_biases), 1 - tolerance, 1 + tolerance)),
      msg = "lapse biases must sum to 1."
    )
  }

  if (lapse_treatment == "no_lapses") {
    lapse_rate <- 0
    lapse_biases <- rep(1 / length(ms), length(ms))
  }

  D <- if (is.list(ms)) length(ms[[1]]) else if (is.matrix(ms)) ncol(ms) else length(ms)
  .assert_that(
    nus[[1]] >= D,
    msg = "Nu must be at least K (number of dimensions of the multivariate Gaussian category)."
  )

  f <- function(x, target_category = 1, logit = FALSE) {
    log_p <- matrix(nrow = length(x), ncol = n.cat)
    for (cat in 1:n.cat) {
      log_p[, cat] <-
        get_NIW_posterior_predictive(
          x,
          ms[[cat]], Ss[[cat]], kappas[[cat]], nus[[cat]],
          Sigma_noise = Sigma_noise, noise_treatment = noise_treatment,
          log = TRUE
        )
    }

    p_target <-
      (1 - lapse_rate) *
        exp(log_p[, target_category] + log(priors[target_category]) -
          log(rowSums(
            t(apply(exp(log_p), 1, function(row) {
              row * priors
            }))
          ))) +
        lapse_rate * lapse_biases[target_category]

    if (logit) qlogis(p_target) else p_target
  }

  return(f)
}

#' Deprecated: get_categorization_function_from_NIW_ideal_adaptor
#'
#' @description `r lifecycle::badge("deprecated")`
#' `get_categorization_function_from_NIW_ideal_adaptor()` was deprecated in MVBeliefUpdatr 0.1.0 and will be
#' removed in 0.2.0. Please use [get_category_posterior_function()] instead.
#'
#' @param model A model object.
#' @param ... Additional arguments.
#' @return A categorization function.
#' @seealso [get_category_posterior_function()]
#' @keywords internal
#' @rdname get_NIW_categorization_function
#' @export
get_categorization_function_from_NIW_ideal_adaptor <- function(model, ...) {
  lifecycle::deprecate_warn(
    "0.1.0",
    "get_categorization_function_from_NIW_ideal_adaptor()",
    with = "get_category_posterior_function()"
  )

  if (S7::S7_inherits(model, MVBU_CognitiveModel)) {
    return(get_category_posterior_function(model, ...))
  }

  suppressWarnings(get_NIW_categorization_function(
    ms = model$m,
    Ss = model$S,
    kappas = model$kappa,
    nus = model$nu,
    priors = model$prior,
    lapse_rate = model$lapse_rate[[1]],
    lapse_biases = model$lapse_bias,
    Sigma_noise = model$Sigma_noise[[1]],
    ...
  ))
}

#' Deprecated: get_categorization_from_NIW_ideal_adaptor
#'
#' @description `r lifecycle::badge("deprecated")`
#' `get_categorization_from_NIW_ideal_adaptor()` was deprecated in MVBeliefUpdatr 0.1.0 and will be
#' removed in 0.2.0. Please use [categorize()] instead.
#'
#' @param x Observations.
#' @param model Model object.
#' @param decision_rule Decision rule.
#' @param noise_treatment Noise treatment.
#' @param lapse_treatment Lapse treatment.
#' @param simplify Logical; whether to simplify output.
#' @param verbose Logical; verbosity.
#' @return Categorization data frame or vector.
#' @seealso [categorize()]
#' @keywords internal
#' @rdname get_categorization_from_model
#' @export
get_categorization_from_NIW_ideal_adaptor <- function(
  x,
  model,
  decision_rule = "sampling",
  noise_treatment = if (decision_rule == "sampling") "sample" else "no_noise",
  lapse_treatment = if (decision_rule == "sampling") "sample" else "marginalize",
  simplify = FALSE,
  verbose = FALSE
) {
  lifecycle::deprecate_warn(
    "0.1.0",
    "get_categorization_from_NIW_ideal_adaptor()",
    with = "categorize()"
  )

  d.response <-
    .legacy_long_posterior(model, x, noise_treatment, lapse_treatment) %>%
    .legacy_apply_decision_rule(decision_rule)

  if (simplify) .legacy_simplify_categorization(d.response, decision_rule) else d.response
}

