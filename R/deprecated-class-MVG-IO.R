get_expected_columns_for_MVG_ideal_observer <- function() append(get_expected_columns_for_MVG(), get_expected_columns_for_model())

#' Deprecated: is.MVG_ideal_observer
#'
#' @description `r lifecycle::badge("deprecated")`
#' `is.MVG_ideal_observer()` was deprecated in MVBeliefUpdatr 0.1.0 and will be removed in 0.2.0.
#' Please use S7 validators and [MVG_IdealObserver] instead.
#'
#' @param x Object to be checked.
#' @param group Name of one or more group variables, each unique combination of which describes an MVG_ideal_observer. (default: NULL)
#' @param category Name of the category variable. (default: "category")
#' @param is.long Is this check assessing whether the ideal observer is in long format (`TRUE`) or wide format (`FALSE`)?
#' (default: `TRUE`)
#' @param with.lapse Does this ideal observer have a lapse rate? (default: `FALSE`)
#' @param with.lapse_bias Does this ideal observer have a lapse bias? (default: `FALSE`)
#' @param verbose Should verbose output be provided? (default: `TRUE`)
#' @param tolerance Probability tolerance.
#'
#' @return A logical.
#'
#' @seealso [MVG_IdealObserver]
#' @keywords internal
#' @export
is.MVG_ideal_observer <- function(x, group = NULL, category = "category", is.long = T, with.lapse = if (with.lapse_bias) T else F, with.lapse_bias = F, verbose = F, tolerance = MVBU_PROB_TOL) {
  lifecycle::deprecate_warn(
    when = "0.1.0",
    what = "is.MVG_ideal_observer()",
    details = "Use S7 validators and MVG_IdealObserver instead."
  )
  name_of_x <- deparse(substitute(x))
  .assert_logical_scalar(with.lapse)
  .assert_logical_scalar(with.lapse_bias)

  if (S7::S7_inherits(x, MVBU_Object)) {
    return(S7::S7_inherits(x, MVG_IdealObserver))
  }

  if (!is.MVBU_model(x, group = group, verbose = verbose, tolerance = tolerance)) {
    return(FALSE)
  }

  # When no groups are specified, infer groups from object.
  if (is.null(group)) {
    group <- setdiff(names(x), get_expected_columns_for_MVG_ideal_observer())
    if (length(group) == 0) group <- NULL else {
      if (verbose) message(paste(name_of_x, "has additional columns beyond those expected:", paste(group, collapse = ", "), "Interpreting those columns as group variables."))
    }
  }


  if (!is.MVG(x, category = category, group = group, verbose = verbose)) {
    if (verbose) message("x does not contain multivariate Gaussian categories.")
    return(FALSE)
  }

  # Only need to test for MVG columns here since is.MVBU_model is called below.
  if (any(get_expected_columns_for_MVG() %nin% names(x))) {
    if (verbose) message(paste("x is missing a required column: ", paste(get_expected_columns_for_MVG, collapse = ",")))
    return(FALSE)
  }

  if (any(!is.factor(get(category, x)))) return(FALSE)

  return(TRUE)
}


