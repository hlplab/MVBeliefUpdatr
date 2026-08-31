#' @include internal-asserts.R
NULL

#' Assert that an object inherits from one of the new S7 classes.
#'
#' @description These helpers validate whether an object inherits from the
#'   S7 classes introduced for category representations, templates, cognitive
#'   models, model distributions, and ideal-adaptor Stanfit objects.
#' @param x Object to test.
#' @param msg Optional error message.
#' @return Invisibly TRUE if the assertion passes.
#' @rdname assert_mvbu_s7_classes
#' @export
assert_MVBU_CategoryRepresentation <- function(x, msg = NULL) {
  .assert_true(
    S7::S7_inherits(x, MVBU_CategoryRepresentation),
    msg = if (is.null(msg)) paste(deparse(substitute(x)), "must inherit from MVBU_CategoryRepresentation") else msg
  )
}

#' @rdname assert_mvbu_s7_classes
#' @export
assert_MVBU_CategoryRepresentationTemplate <- function(x, msg = NULL) {
  .assert_true(
    S7::S7_inherits(x, MVBU_CategoryRepresentationTemplate),
    msg = if (is.null(msg)) paste(deparse(substitute(x)), "must inherit from MVBU_CategoryRepresentationTemplate") else msg
  )
}

#' @rdname assert_mvbu_s7_classes
#' @export
assert_MVBU_CognitiveModel <- function(x, msg = NULL) {
  .assert_true(
    S7::S7_inherits(x, MVBU_CognitiveModel),
    msg = if (is.null(msg)) paste(deparse(substitute(x)), "must inherit from MVBU_CognitiveModel") else msg
  )
}

#' @rdname assert_mvbu_s7_classes
#' @export
assert_IdealAdaptorStanfit <- function(x, msg = NULL) {
  .assert_true(
    S7::S7_inherits(x, IdealAdaptorStanfit),
    msg = if (is.null(msg)) paste(deparse(substitute(x)), "must inherit from IdealAdaptorStanfit") else msg
  )
}

#' @rdname assert_mvbu_s7_classes
#' @export
assert_IdealAdaptorStaninput <- function(x, msg = NULL) {
  .assert_true(S7::S7_inherits(x, IdealAdaptorStaninput), msg = if (is.null(msg)) "Expected an IdealAdaptorStaninput object." else msg)
}
