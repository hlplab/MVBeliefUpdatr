#' @include internal-asserts.R
NULL

#' @name assert_MVG_ideal_observer
#' @title Deprecated assertion helper for MVG ideal observer objects
#' @description Deprecated. This compatibility helper is deprecated and will be removed in a future release.
#' @param x Object to test.
#' @param category Name of the category field to inspect.
#' @param verbose Logical. Whether to emit additional diagnostics while testing the object.
#' @return Invisibly TRUE if the assertion passes.
#' @keywords internal
#' @export
assert_MVG_ideal_observer = function(x, category = "category", verbose = F) {
  .assert_that(is.MVG_ideal_observer(x, category = category, verbose = verbose),
               msg = paste(deparse(substitute(x)), "must be an MVG_ideal_observer object."))
}

#' @name assert_NIW_belief
#' @title Deprecated assertion helper for NIW belief objects
#' @description Deprecated. Use \code{\link{assert_IdealAdaptorStanfit}} for Stanfit objects, \code{\link{assert_staninput}} for Staninput objects, or \code{\link{assert_stanfit_input}} for Stanfit input containers instead. This compatibility helper is deprecated and will be removed in a future release.
#' @keywords internal
#' @param x Object to test.
#' @param category Name of the category field to inspect.
#' @param verbose Logical. Whether to emit additional diagnostics while testing the object.
#' @param strict Logical. Whether to require the object to be an NIW belief strictly, rather than allowing an ideal adaptor.
#' @return Invisibly TRUE if the assertion passes.
#' @export
assert_NIW_belief = function(x, category = "category", verbose = F, strict = F) {
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

#' @name assert_NIW_ideal_adaptor
#' @title Deprecated assertion helper for NIW ideal adaptor objects
#' @description Deprecated. Use \code{\link{assert_IdealAdaptorStanfit}} for Stanfit objects, \code{\link{assert_staninput}} for Staninput objects, or \code{\link{assert_stanfit_input}} for Stanfit input containers instead. This compatibility helper is deprecated and will be removed in a future release.
#' @keywords internal
#' @param x Object to test.
#' @param category Name of the category field to inspect.
#' @param verbose Logical. Whether to emit additional diagnostics while testing the object.
#' @return Invisibly TRUE if the assertion passes.
#' @export
assert_NIW_ideal_adaptor = function(x, category = "category", verbose = F) {
  .assert_that(is.NIW_ideal_adaptor(x, category = category, verbose = verbose),
               msg = paste(deparse(substitute(x)), "must be an NIW_ideal_adaptor object."))
}

