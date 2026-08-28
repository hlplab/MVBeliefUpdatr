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
#'   For details, see \code{\link{get_NIW_posterior_predictive}}. Note though that "sample" would likely result in rather
#'   counter-intuitive behavior of the categorization function and is thus not recommended.
#' @param lapse_treatment Should the consequences of attentional lapses be included in the categorization function
#'   ("marginalize") or not ("no_lapses")? (default: "marginalize")
#'
#' @return A function that takes as input cue values and returns posterior probabilities of the first category,
#'   based on the posterior predictive of the cues given the (IBBU-derived parameters for the) categories' m, S,
#'   kappa, nu, and prior, as well as the lapse rate. The function will accept the following arguments:
#'   \itemize{
#'     \item{`target_category`:}{The index of the category for which categorization should be shown. (default: `1`)}
#'     \item{`logit`:}{Should the function return log-odds (TRUE) or probabilities (FALSE)? (default: FALSE)}
#'   }
#'
#' @seealso TBD
#'
#' @rdname get_NIW_categorization_function
#' @importFrom purrr map_lgl
#' @export
get_NIW_categorization_function <- function(
    ms, Ss, kappas, nus,
    priors = rep(1 / length(ms), length(ms)),
    lapse_rate = 0,
    lapse_biases = rep(1 / length(ms), length(ms)),
    Sigma_noise = NULL,
    noise_treatment = infer_default_noise_treatment(Sigma_noise),
    lapse_treatment = if (lapse_rate > 0) "marginalize" else "no_lapses"
) {
  tolerance = MVBU_PROB_TOL
  .assert_that(.is_equal(length(ms), length(Ss)),
              .is_equal(length(ms), length(priors)),
              .is_equal(length(ms), length(kappas)),
              .is_equal(length(ms), length(nus)),
              msg = "The number of ms, Ss, kappas, nus, and priors must be identical.")
  n.cat = length(ms)

  .assert_that(all(between(priors, 0, 1), between(sum(priors), 1 - tolerance, 1 + tolerance)),
              msg = "priors must sum to 1.")
  .assert_that(.is_non_NA_scalar_double(lapse_rate),
              msg = "lapse_rate must be a scalar.")
  .assert_that(between(lapse_rate, 0, 1))
  if (any(is.null(lapse_biases),
          all(is.null(lapse_biases)),
          all(map_lgl(lapse_biases, is.null)))) {
    lapse_biases <- 1 / n.cat
  } else {
    .assert_that(all(between(lapse_biases, 0, 1), between(sum(lapse_biases), 1 - tolerance, 1 + tolerance)),
                msg = "lapse biases must sum to 1.")
  }

  if (lapse_treatment == "no_lapses") {
    lapse_rate = 0
    lapse_biases = rep(1 / length(ms), length(ms))
  }

  # Get dimensions of multivariate category
  D = get_D(ms)
  .assert_that(
    nus[[1]] >= D,
    msg = "Nu must be at least K (number of dimensions of the multivariate Gaussian category).")

  f <- function(x, target_category = 1, logit = F) {
    log_p <- matrix(nrow = length(x), ncol = n.cat)
    for (cat in 1:n.cat) {
      log_p[, cat] <-
        get_NIW_posterior_predictive(
          x,
          ms[[cat]], Ss[[cat]], kappas[[cat]], nus[[cat]],
          Sigma_noise = Sigma_noise, noise_treatment = noise_treatment,
          log = T)
    }

    p_target <-
      (1 - lapse_rate) *
      exp(log_p[,target_category] + log(priors[target_category]) -
            log(rowSums(
              # multiply each row of the exponentiated log_p matrix (= the log densities of each category for each observation)
              # with the vector of priors (must use this since simply multiplying by the vector will try to proceed element-wise
              # through the matrix, through each column)
              t(apply(exp(log_p), 1, function(row) { row * priors }))))) +
      lapse_rate * lapse_biases[target_category]

    if (logit) return(qlogis(p_target)) else return(p_target)
  }

  return(f)
}

# Deprecated after S7-migration

#' Legacy wrapper for get_category_posterior_function.
#'
#' @description Deprecated. Use \code{\link{get_category_posterior_function}} instead.
#' @rdname get_NIW_categorization_function
#' @export
#' @description Deprecated. Use get_category_posterior_function() instead.
#' @keywords internal
get_categorization_function_from_NIW_ideal_adaptor <- function(model, ...) {
  lifecycle::deprecate_warn(
    when = "0.0.3",
    what = "get_categorization_function_from_NIW_ideal_adaptor()",
    with = "get_category_posterior_function()"
  )
  # Could be used later in a function that checks internal consistency of model
  # if (nunique(model$lapse_rate) > 1) stop2("Model has more than one unique lapse_rate.")
  #
  # .is_identical_to_first <- function(x, first_matrix) {
  #   identical(x, first_matrix)
  # }
  # if (nunique(model$lapse_rate) > 1) stop2("Model has more than one unique Sigma_noise.")

  get_NIW_categorization_function(
    ms = model$m,
    Ss = model$S,
    kappas = model$kappa,
    nus = model$nu,
    priors = model$prior,
    lapse_rate = model$lapse_rate[[1]],
    lapse_biases = model$lapse_bias,
    Sigma_noise = model$Sigma_noise[[1]],
    ...
  )
}

#' Legacy wrapper for categorize.
#'
#' @description Deprecated. Use \code{\link{categorize}} instead.
#' @rdname get_categorization_from_model
#' @export
#' @description Deprecated. Use categorize() instead.
#' @keywords internal
get_categorization_from_NIW_ideal_adaptor <- function(
  x,
  model,
  decision_rule = "sampling",
  noise_treatment = if (decision_rule == "sampling") "sample" else "no_noise",
  lapse_treatment = if (decision_rule == "sampling") "sample" else "marginalize",
  simplify = F,
  verbose = F
) {
  lifecycle::deprecate_warn(
    when = "0.0.3",
    what = "get_categorization_from_NIW_ideal_adaptor()",
    with = "categorize()"
  )

  d.response <-
    .legacy_long_posterior(model, x, noise_treatment, lapse_treatment) %>%
    .legacy_apply_decision_rule(decision_rule)

  if (simplify) .legacy_simplify_categorization(d.response, decision_rule) else d.response
}
