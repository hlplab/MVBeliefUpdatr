#' @include S7-core-niw-classes.R
#' @importFrom lifecycle deprecate_warn
#' @importFrom S7 S7_inherits
NULL

#' Deprecated: An S4 class for legacy NIW ideal adaptor objects
#'
#' @name NIW_ideal_adaptor-class
#' @aliases NIW_ideal_adaptor
#' @keywords internal
#' @noRd
NIW_ideal_adaptor <-
  setClass(
    "NIW_ideal_adaptor",
    contains = "tbl_df",
    package = "MVBeliefUpdatr"
  )

# Call class constructor function
NIW_ideal_adaptor

#' @keywords internal
#' @noRd
get_expected_columns_for_NIW_ideal_adaptor <- function() {
  c(get_expected_columns_for_NIW_belief(), get_expected_columns_for_model())
}

# deprecated ------------------------------------------------------------------

#' Deprecated: is.NIW_ideal_adaptor
#'
#' @description `r lifecycle::badge("deprecated")`
#' `is.NIW_ideal_adaptor()` was deprecated in MVBeliefUpdatr 0.1.0 and will be removed in 0.2.0.
#' Please use S7 validators and [NIW_IdealAdaptor] instead.
#'
#' @param x Object to check.
#' @param group Name of one or more group variables. (default: NULL)
#' @param category Name of the category variable. (default: "category")
#' @param is.long Is this check assessing whether the ideal adaptor is in long format? (default: `TRUE`)
#' @param with.prior Does this ideal adaptor have a prior? (default: `TRUE`)
#' @param with.lapse Does this ideal adaptor have a lapse rate? (default: `FALSE`)
#' @param with.lapse_bias Does this ideal adaptor have a lapse bias? (default: `FALSE`)
#' @param verbose Should verbose output be provided? (default: `FALSE`)
#' @param tolerance Tolerance for sum-to-one probability checks. (default: MVBU_PROB_TOL)
#' @param ... Additional arguments.
#' @return Logical indicating whether `x` is a valid NIW ideal adaptor.
#' @seealso [NIW_IdealAdaptor]
#' @keywords internal
#' @export
is.NIW_ideal_adaptor <- function(x, group = NULL, category = "category", is.long = T, with.prior = T, with.lapse = if (with.lapse_bias) T else F, with.lapse_bias = F, verbose = F, tolerance = MVBU_PROB_TOL, ...) {
  lifecycle::deprecate_warn(
    when = "0.1.0",
    what = "is.NIW_ideal_adaptor()",
    details = "Use S7::S7_inherits(x, NIW_IdealAdaptor) instead."
  )
  name_of_x <- deparse(substitute(x))
  .assert_logical_scalar(with.lapse)
  .assert_logical_scalar(with.lapse_bias)

  if (S7::S7_inherits(x, MVBU_Object)) {
    return(S7::S7_inherits(x, NIW_IdealAdaptor))
  }

  if (!is.MVBU_model(x, group = group, verbose = verbose, tolerance = tolerance)) {
    return(FALSE)
  }

  # When no groups are specified, infer groups from object.
  if (is.null(group)) {
    group <- setdiff(names(x), get_expected_columns_for_NIW_ideal_adaptor())
    if (length(group) == 0) group <- NULL else {
      if (verbose) message(paste(name_of_x, "has additional columns beyond those expected:", paste(group, collapse = ", "), "Interpreting those columns as group variables."))
    }
  }

  if (!is.null(group)) {
    if (verbose) message("Checking whether ", name_of_x, " is an NIW_ideal_adaptor within each unique combination of group values.")
    x %<>% group_by(!!! syms(group))
  }

  if (!is.NIW_belief(x, group = group)) {
    if (verbose) message(paste(deparse(substitute(x)), "does not contain NIW beliefs."))
    return(FALSE)
  }

  if (
    any(
      !with.prior | "prior" %nin% names(x),
      with.lapse & "lapse_rate" %nin% names(x),
      with.lapse_bias & "lapse_bias" %nin% names(x)
    )
  ) {
    if (verbose) message(paste(name_of_x, " is missing prior, lapse rate, or lapse bias."))
    return(FALSE)
  }

  return(TRUE)
}
