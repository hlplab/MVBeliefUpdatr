get_expected_columns_for_exemplar_model <- function() append(get_expected_columns_for_exemplars(), get_expected_columns_for_model())

#' Deprecated: is.exemplar_model
#'
#' @description `r lifecycle::badge("deprecated")`
#' `is.exemplar_model()` was deprecated in MVBeliefUpdatr 0.1.0 and will be removed in 0.2.0.
#' Please use S7 validators and [Exemplar_Model] instead.
#'
#' @param x Object to be checked.
#' @param group Name of one or more group variables, each unique combination of which describes an exemplar model. (default: NULL)
#' @param verbose Should verbose output be provided? (default: `TRUE`)
#' @param tolerance Probability tolerance.
#'
#' @return A logical.
#'
#' @seealso [Exemplar_Model]
#' @keywords internal
#' @export
is.exemplar_model <- function(x, group = NULL, verbose = F, tolerance = MVBU_PROB_TOL) {
  lifecycle::deprecate_warn(
    when = "0.1.0",
    what = "is.exemplar_model()",
    details = "Use S7 validators and Exemplar_Model instead."
  )
  name_of_x <- deparse(substitute(x))

  if (!is.MVBU_model(x, group = group, verbose = verbose, tolerance = tolerance)) {
    return(FALSE)
  }

  # When no groups are specified, infer groups from object.
  if (is.null(group)) {
    group <- setdiff(names(x), get_expected_columns_for_exemplar_model())
    if (length(group) == 0) group <- NULL else {
      if (verbose) message(name_of_x, "has additional columns beyond those expected:", paste(group, collapse = ", "), "Interpreting those columns as group variables.")
    }
  }

  if (!is.exemplars(x, group = group, verbose = verbose)) {
    if (verbose) message("x does not contain exemplars.")
    return(FALSE)
  }

  # Only need to test for exemplars columns here since is.MVBU_model is called below.
  if (any(get_expected_columns_for_exemplars() %nin% names(x))) {
    if (verbose) message("x is missing a required column: ", paste(get_expected_columns_for_exemplars, collapse = ","))
    return(FALSE)
  }

  if (any(!is.factor(x$category))) {
    if (verbose) message(paste("category must be a factor."))
    return(FALSE)
  }

  return(TRUE)
}


