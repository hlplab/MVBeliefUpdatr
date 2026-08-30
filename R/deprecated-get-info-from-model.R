#' @include S7-core-classes.R
#' @importFrom lifecycle deprecate_warn
#' @importFrom tibble tibble
NULL

# Internal legacy helpers

.legacy_observation_matrix <- function(x, d, arg_name = "x") {
  if (is.list(x) && !is.data.frame(x)) {
    if (length(x) == 0L) .stop(arg_name, " must contain at least one observation.")
    return(do.call(rbind, lapply(x, function(v) .as_observation_matrix(v, d = d, arg_name = arg_name))))
  }
  .as_observation_matrix(x, d = d, arg_name = arg_name)
}

.legacy_observation_list <- function(x_mat) {
  lapply(seq_len(nrow(x_mat)), function(i) unname(x_mat[i, ]))
}

.legacy_long_posterior <- function(model, x, noise_treatment, lapse_treatment) {
  d <- length(get_cue_labels(model))
  x_mat <- .legacy_observation_matrix(x, d)
  x_list <- .legacy_observation_list(x_mat)

  pf <- get_category_posterior_function(model, noise_treatment, lapse_treatment)
  post <- pf(x_mat, categories = NULL)
  cats <- colnames(post)

  tibble::tibble(
    observationID = rep(seq_len(nrow(post)), each = length(cats)),
    x = rep(x_list, each = length(cats)),
    category = rep(cats, times = nrow(post)),
    posterior_probability = as.vector(t(post))
  )
}

.legacy_apply_decision_rule <- function(d.post, decision_rule) {
  n_cats <- length(unique(d.post$category))
  probs <- matrix(d.post$posterior_probability, ncol = n_cats, byrow = TRUE)

  response <- if (identical(decision_rule, "proportional")) {
    probs
  } else if (identical(decision_rule, "criterion")) {
    t(apply(probs, 1, function(p) as.numeric(seq_along(p) == which.max(p))))
  } else if (identical(decision_rule, "sampling")) {
    t(apply(probs, 1, function(p) as.numeric(stats::rmultinom(1, 1, p))))
  } else {
    .stop("Decision rule must be one of: 'criterion', 'proportional', or 'sampling'.")
  }

  d.post$response <- as.vector(t(matrix(response, ncol = n_cats)))
  d.post[, c("observationID", "x", "category", "response")]
}

.legacy_simplify_categorization <- function(d.response, decision_rule) {
  .assert_that(decision_rule %in% c("criterion", "sampling"),
    msg = "For simplify = T, decision rule must be either criterion or sampling."
  )
  winners <- d.response[d.response$response == 1, c("observationID", "category")]
  winners <- winners[order(winners$observationID), ]
  winners$category
}

# Deprecated functions

#' Get cue dimensionality from likelihood or model (deprecated)
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#' `get_cue_dimensionality_from_model()` is deprecated. Use `length(get_cue_labels(x))` instead.
#'
#' @param x A likelihood or model object.
#' @param indices Optional indices (deprecated and ignored).
#'
#' @return A numeric integer.
#' @keywords internal
#' @export
get_cue_dimensionality_from_model <- function(x, indices = NULL) {
  lifecycle::deprecate_warn(
    "0.1.0",
    "get_cue_dimensionality_from_model()",
    details = "Use length(get_cue_labels(x)) instead."
  )
  length(get_cue_labels(x))
}

#' Get noise covariance matrix from model (deprecated)
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#' `get_perceptual_noise_from_model()` is deprecated. Use `get_noise()` instead.
#'
#' @param model A model object.
#' @return Perceptual noise covariance matrix.
#' @keywords internal
#' @export
get_perceptual_noise_from_model <- function(model) {
  lifecycle::deprecate_warn(
    "0.1.0",
    "get_perceptual_noise_from_model()",
    "get_noise()"
  )
  get_noise(model)
}

#' Legacy wrapper for posterior.
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#' `get_posterior_from_model()` is deprecated. Use `posterior()` instead.
#'
#' @param model A model object.
#' @param ... Arguments passed to \code{\link{posterior}}.
#'
#' @rdname get_posterior_from_model
#' @export
#' @keywords internal
get_posterior_from_model <- function(model, ...) {
  lifecycle::deprecate_warn(
    when = "0.1.0",
    what = "get_posterior_from_model()",
    with = "posterior()",
    always = TRUE
  )
  dots <- list(...)
  if (!is.null(dots$x) && is.null(dots$new_data)) {
    dots$new_data <- dots$x
    dots$x <- NULL
  }
  do.call(posterior, c(list(model), dots))
}

#' Legacy wrapper for categorize.
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#' `get_categorization_from_model()` is deprecated. Use `categorize()` instead.
#'
#' @param model A model object.
#' @param decision_rule Decision rule to use.
#' @param ... Arguments passed to \code{\link{categorize}}.
#'
#' @rdname get_categorization_from_model
#' @export
#' @keywords internal
get_categorization_from_model <- function(model, decision_rule = "sampling", ...) {
  lifecycle::deprecate_warn(
    when = "0.1.0",
    what = "get_categorization_from_model()",
    with = "categorize()",
    always = TRUE
  )
  dots <- list(...)
  if (!is.null(dots$x) && is.null(dots$new_data)) {
    dots$new_data <- dots$x
    dots$x <- NULL
  }
  do.call(categorize, c(list(model), dots, list(decision_rule = decision_rule)))
}

#' Get cue labels from likelihood or model (deprecated)
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#' `get_cue_labels_from_model()` is deprecated. Use `get_cue_labels()` instead.
#'
#' @param x A model, template, or representation.
#' @param indices Optional integer indices of cue labels to return.
#' @return Character vector of cue labels.
#' @keywords internal
#' @export
get_cue_labels_from_model <- function(x, indices = NULL) {
  lifecycle::deprecate_warn(
    "0.1.0",
    "get_cue_labels_from_model()",
    "get_cue_labels()"
  )
  get_cue_labels(x, indices = indices)
}

#' Get category labels from likelihood or model (deprecated)
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#' `get_category_labels_from_model()` is deprecated. Use `get_category_labels()` instead.
#'
#' @param x A model, template, or representation.
#' @param indices Optional integer indices of category labels to return.
#' @return Character vector of category labels.
#' @keywords internal
#' @export
get_category_labels_from_model <- function(x, indices = NULL) {
  lifecycle::deprecate_warn(
    "0.1.0",
    "get_category_labels_from_model()",
    "get_category_labels()"
  )
  get_category_labels(x, indices = indices)
}

#' Get number of categories from likelihood or model (deprecated)
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#' `get_nlevels_of_category_labels_from_model()` is deprecated. Use `length(get_category_labels())` instead.
#'
#' @param x A model, template, or representation.
#' @return Number of categories as an integer.
#' @keywords internal
#' @export
get_nlevels_of_category_labels_from_model <- function(x) {
  lifecycle::deprecate_warn(
    "0.1.0",
    "get_nlevels_of_category_labels_from_model()",
    details = "Use length(get_category_labels(x)) instead."
  )
  length(get_category_labels(x))
}

#' Get category priors from model (deprecated)
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#' `get_priors_from_model()` is deprecated. Use `get_category_prior()` instead.
#'
#' @param model A model object.
#' @param categories Optional vector of category labels to subset.
#' @return Numeric vector of prior probabilities.
#' @keywords internal
#' @export
get_priors_from_model <- function(model, categories = NULL) {
  lifecycle::deprecate_warn(
    "0.1.0",
    "get_priors_from_model()",
    "get_category_prior()"
  )
  prior <- get_category_prior(model)
  if (!is.null(categories)) {
    prior <- prior[match(as.character(categories), names(prior))]
  }
  as.numeric(prior)
}

#' Get lapse rate from model (deprecated)
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#' `get_lapse_rate_from_model()` is deprecated. Use `get_lapse_rate()` instead.
#'
#' @param model A model object.
#' @return Numeric lapse rate.
#' @keywords internal
#' @export
get_lapse_rate_from_model <- function(model) {
  lifecycle::deprecate_warn(
    "0.1.0",
    "get_lapse_rate_from_model()",
    "get_lapse_rate()"
  )
  as.numeric(get_lapse_rate(model))
}

#' Get lapse biases from model (deprecated)
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#' `get_lapse_biases_from_model()` is deprecated. Use `get_lapse_bias()` instead.
#'
#' @param model A model object.
#' @param categories Optional vector of category labels to subset.
#' @return Numeric vector of lapse biases.
#' @keywords internal
#' @export
get_lapse_biases_from_model <- function(model, categories = NULL) {
  lifecycle::deprecate_warn(
    "0.1.0",
    "get_lapse_biases_from_model()",
    "get_lapse_bias()"
  )
  bias <- get_lapse_bias(model)
  if (!is.null(categories)) {
    bias <- bias[match(as.character(categories), names(bias))]
  }
  as.numeric(bias)
}
