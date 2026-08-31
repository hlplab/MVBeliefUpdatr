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
#' @description
#' `r lifecycle::badge("deprecated")`
#' `is.NIW_ideal_adaptor()` is deprecated. Use
#' `S7::S7_inherits(x, NIW_IdealAdaptor)` instead.
#'
#' @param x Object to check.
#' @param ... Additional arguments (ignored; for compatibility).
#' @return Logical indicating whether `x` inherits from [NIW_IdealAdaptor].
#' @rdname deprecated-functions
#' @export
is.NIW_ideal_adaptor <- function(x, ...) {
  lifecycle::deprecate_warn(
    when = "0.2.0",
    what = "is.NIW_ideal_adaptor()",
    details = "Use S7::S7_inherits(x, NIW_IdealAdaptor) instead."
  )
  S7::S7_inherits(x, NIW_IdealAdaptor)
}
