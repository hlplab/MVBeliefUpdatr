NULL

#' @name get_class
#' @title Deprecated class-dispatch helper
#' @description Deprecated. Use the S7-based validators and constructors instead.
#' @keywords internal
#' @export
get_class <- function(x) {
  lifecycle::deprecate_warn(
    when = "0.0.3",
    what = "get_class()",
    details = "Use the S7-based validators and constructors."
  )
  if (is.NIW_ideal_adaptor(x)) {
    return("NIW_ideal_adaptor")
  }
  if (is.NIW_belief(x)) {
    return("NIW_belief")
  }
  if (is.MVG_ideal_observer(x)) {
    return("MVG_ideal_observer")
  }
  if (is.MVG(x)) {
    return("MVG")
  }
  if (is.exemplar_model(x)) {
    return("exemplar_model")
  }
  if (is.exemplars(x)) {
    return("exemplars")
  }
  if (is.ideal_adaptor_stanfit(x)) {
    return("ideal_adaptor_stanfit")
  }

  "Unrecognized class"
}

#' @name get_expected_columns_for_model
#' @title Deprecated expected-column helper
#' @description Deprecated. Use the S7-based predicates and constructors instead.
#' @keywords internal
get_expected_columns_for_model <- function() {
  c("prior", "lapse_rate", "lapse_bias", "Sigma_noise")
}

.split_by_group <- function(x, group) {
  if (is.null(group) || length(group) == 0L) {
    return(list(seq_len(nrow(x))))
  }

  if (!all(group %in% names(x))) {
    return(list(seq_len(nrow(x))))
  }

  group_key <- interaction(x[group], drop = TRUE)
  split(seq_len(nrow(x)), group_key)
}

#' Is this an MVBeliefUpdatr representation?
#'
#' Check whether \code{x} is recognized as an MVBeliefUpdatr category representation.
#'
#' @param x Object to be checked.
#' @param group Name of one or more group variables, each unique combination of which describes a model. (default: NULL)
#' @param verbose Should verbose output be provided? (default: `TRUE`)
#'
#' @return A logical.
#'
#' @seealso TBD
#' @description Deprecated. Use the S7-based predicates and constructors instead.
#' @keywords internal
#' @export
is.MVBU_representation <- function(x, group = NULL, verbose = F, tolerance = MVBU_PROB_TOL) {
  lifecycle::deprecate_warn(
    when = "0.0.3",
    what = "is.MVBU_representation()",
    details = "Use the S7-based predicates and constructors."
  )
  name_of_x <- deparse(substitute(x))

  if (!is.data.frame(x)) {
    if (verbose) message("Object is not a data frame-like object.")
    return(FALSE)
  }

  if (!is.null(group)) {
    if (verbose) message("Checking whether ", name_of_x, " is a model within each unique combination of group values.")
    if (!all(group %in% names(x))) {
      return(FALSE)
    }
  }

  rows <- .split_by_group(x, group)
  for (idx in rows) {
    subset <- x[idx, , drop = FALSE]
    if (!all(vapply(seq_len(nrow(subset)), function(i) {
      is.exemplars(subset[i, , drop = FALSE], verbose = FALSE) ||
        is.MVG(subset[i, , drop = FALSE], verbose = FALSE) ||
        is.NIW_belief(subset[i, , drop = FALSE], verbose = FALSE)
    }, logical(1)))) {
      return(FALSE)
    }
  }

  TRUE
}

#' Is this an MVBeliefUpdatr model?
#'
#' Check whether \code{x} is a model with lapse rates, biases, priors, and perceptual noise.
#'
#' @param x Object to be checked.
#' @param group Name of one or more group variables, each unique combination of which describes a model. (default: NULL)
#' @param verbose Should verbose output be provided? (default: `TRUE`)
#'
#' @return A logical.
#'
#' @seealso TBD
#' @description Deprecated. Use the S7-based predicates and constructors instead.
#' @keywords internal
#' @export
is.MVBU_model <- function(x, group = NULL, verbose = F, tolerance = MVBU_PROB_TOL) {
  lifecycle::deprecate_warn(
    when = "0.0.3",
    what = "is.MVBU_model()",
    details = "Use the S7-based validators and constructors."
  )
  name_of_x <- deparse(substitute(x))

  if (!is.data.frame(x)) {
    if (verbose) message("Object is not a data frame-like object. All MVBeliefUpdatr models are stored in data frames.")
    return(FALSE)
  }

  if (verbose) message("Checking whether ", name_of_x, " has all the column names required for a model.")
  if (!all(c("prior", "lapse_rate", "lapse_bias", "Sigma_noise") %in% names(x))) {
    return(FALSE)
  }


  sigma_noise <- x[["Sigma_noise"]]
  if (!all(vapply(sigma_noise, function(value) is.null(value) || is.matrix(value), logical(1)))) {
    if (verbose) message("If not NULL, Sigma_noise must be a matrix.")
    return(FALSE)
  }

  if (length(sigma_noise) > 0L && !all(vapply(sigma_noise, function(value) identical(value, sigma_noise[[1]]), logical(1)))) {
    if (verbose) message(paste("Noise covariance matrix Sigma_noise in", name_of_x, "is not constant across categories."))
    return(FALSE)
  }

  if (!is.null(sigma_noise[[1]])) {
    sigma_noise_value <- sigma_noise[[1]]
    .assert_that(all(dim(sigma_noise_value) == rep(get_cue_dimensionality_from_model(x), 2)),
                 msg = paste("If not NULL, Sigma_noise must match the dimensionality of other parameters in the model (here: a",
                             get_cue_dimensionality_from_model(x), "x", get_cue_dimensionality_from_model(x), " matrix)."))
    .assert_that(!is.null(dimnames(sigma_noise_value)),
                 msg = "If not NULL, Sigma_noise must have non-NULL dimnames.")
    cue_labels <- get_cue_labels_from_model(x)
    dim_names <- dimnames(sigma_noise_value)
    .assert_that(!is.null(dim_names),
                 msg = "If not NULL, the dimnames of Sigma_noise must match the cue names used in the model.")
    .assert_that(length(dim_names) == 2L && all(dim_names[[1]] == cue_labels) && all(dim_names[[2]] == cue_labels),
                 msg = "If not NULL, the dimnames of Sigma_noise must match the cue names used in the model.")
  }

  if (any(x$prior < 0 | x$prior > 1)) {
    if (verbose) message(paste("Prior probabilities in", name_of_x, "are not all between 0 and 1: ", paste(x$prior, collapse = ",")))
    return(FALSE)
  }

  if (abs(sum(x$prior) - 1) > tolerance) {
    if (verbose) message(paste("Prior probabilities in", name_of_x, "do not add up to 1: ", sum(x$prior)))
    return(FALSE)
  }

  if (any(x$lapse_rate < 0 | x$lapse_rate > 1)) {
    if (verbose) message(paste("Lapse rates in", name_of_x, "are not all between 0 and 1: ", paste(x$lapse_rate, collapse = ",")))
    return(FALSE)
  }

  if (length(x$lapse_rate) > 0L && !all(vapply(x$lapse_rate, function(value) identical(value, x$lapse_rate[[1]]), logical(1)))) {
    if (verbose) {
      message(
        paste(
          "Lapse rates in",
          name_of_x,
          "are not constant across categories: ",
          paste(x$lapse_rate, collapse = ", ")))
    }
    return(FALSE)
  }

  if (any(x$lapse_bias < 0 | x$lapse_bias > 1)) {
    if (verbose) message(paste("Lapse bias probabilities in", name_of_x, "are not all between 0 and 1: ", paste(x$lapse_bias, collapse = ",")))
    return(FALSE)
  }

  if (abs(sum(x$lapse_bias) - 1) > tolerance) {
    if (verbose) message(paste("Lapse bias probabilities in", name_of_x, "do not add up to 1: ", sum(x$lapse_bias)))
    return(FALSE)
  }

  TRUE
}

#' Print MVBeliefUpdatr model
#'
#' Specifies reasonable defaults for the parameters to be summarized for the MVBeliefUpdatr_model object.
#'
#' @param x An \code{\link{MVBeliefUpdatr_model}} object.
#'
#' @description Deprecated. Use the S7-based predicates and constructors instead.
#' @keywords internal
#' @export
print.MVBU_model <- function(x, ...) {
  lifecycle::deprecate_warn(
    when = "0.0.3",
    what = "print.MVBU_model()",
    details = "Use the S7-based validators and constructors."
  )
  .assert_that(is.MVBU_model(x), msg = "Expected an MVBU model.")

  if (get_cue_dimensionality_from_model(x) == 1L) {
    x[] <- lapply(x, function(value) {
      if (is.numeric(value)) unlist(value) else value
    })
  }

  print(x)
}
