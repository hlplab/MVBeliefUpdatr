#' @include internal-asserts.R
NULL

#' Deprecated: assert_MVG_ideal_observer
#'
#' @description `r lifecycle::badge("deprecated")`
#' `assert_MVG_ideal_observer()` was deprecated in MVBeliefUpdatr 0.1.0 and will be
#' removed in 0.2.0. Please use S7 assertions or [MVG_IdealObserver] instead.
#' @param x Object to test.
#' @param category Name of the category field to inspect.
#' @param verbose Logical. Whether to emit additional diagnostics while testing the object.
#' @return Invisibly TRUE if the assertion passes.
#' @seealso [MVG_IdealObserver]
#' @keywords internal
#' @export
assert_MVG_ideal_observer <- function(x, category = "category", verbose = F) {
  lifecycle::deprecate_warn("0.1.0", "assert_MVG_ideal_observer()")
  .assert_that(is.MVG_ideal_observer(x, category = category, verbose = verbose),
               msg = paste(deparse(substitute(x)), "must be an MVG_ideal_observer object."))
}

#' Deprecated: assert_NIW_belief
#'
#' @description `r lifecycle::badge("deprecated")`
#' `assert_NIW_belief()` was deprecated in MVBeliefUpdatr 0.1.0 and will be
#' removed in 0.2.0. Please use [NIW_IdealAdaptor] or S7 assertions instead.
#' @keywords internal
#' @param x Object to test.
#' @param category Name of the category field to inspect.
#' @param verbose Logical. Whether to emit additional diagnostics while testing the object.
#' @param strict Logical. Whether to require the object to be an NIW belief strictly, rather than allowing an ideal adaptor.
#' @return Invisibly TRUE if the assertion passes.
#' @seealso [NIW_IdealAdaptor]
#' @export
assert_NIW_belief <- function(x, category = "category", verbose = F, strict = F) {
  lifecycle::deprecate_warn("0.1.0", "assert_NIW_belief()")
  if (strict) {
    .assert_that(is.NIW_belief(x, category = category, verbose = verbose),
                 msg = paste(deparse(substitute(x)), "must be an NIW_belief object."))
  } else {
    .assert_that(
      any(
        is.NIW_belief(x, category = category, verbose = verbose),
        is.NIW_ideal_adaptor(x, category = category, verbose = verbose)),
      msg = paste(deparse(substitute(x)), "must be an NIW_belief or NIW_ideal_adaptor object."))
  }
}

#' Deprecated: assert_NIW_ideal_adaptor
#'
#' @description `r lifecycle::badge("deprecated")`
#' `assert_NIW_ideal_adaptor()` was deprecated in MVBeliefUpdatr 0.1.0 and will be
#' removed in 0.2.0. Please use [NIW_IdealAdaptor] or S7 assertions instead.
#' @keywords internal
#' @param x Object to test.
#' @param category Name of the category field to inspect.
#' @param verbose Logical. Whether to emit additional diagnostics while testing the object.
#' @return Invisibly TRUE if the assertion passes.
#' @seealso [NIW_IdealAdaptor]
#' @export
assert_NIW_ideal_adaptor <- function(x, category = "category", verbose = F) {
  lifecycle::deprecate_warn("0.1.0", "assert_NIW_ideal_adaptor()")
  .assert_that(is.NIW_ideal_adaptor(x, category = category, verbose = verbose),
               msg = paste(deparse(substitute(x)), "must be an NIW_ideal_adaptor object."))
}


