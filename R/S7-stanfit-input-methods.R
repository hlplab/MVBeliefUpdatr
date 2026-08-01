#' @include S7-core-classes.R
#' @include S7-generics.R
#' @include S7-transform-information.R
#' @include S7-staninput.R
#' @include S7-stanfit-input.R
#' @include S7-stanfit.R
#' @include S7-stanfit-methods.R
NULL

S7::method(get_staninput, IdealAdaptorStanfitInput) <- function(x) {
  x@staninput
}

S7::method(set_staninput, list(S7::class_any, S7::class_any)) <- function(x, staninput) {
  stop("x must be an IdealAdaptorStanfit or IdealAdaptorStanfitInput object.", call. = FALSE)
}

S7::method(set_staninput, list(MVBU_Stanfit, S7::class_any)) <- function(x, staninput) {
  if (!S7::S7_inherits(staninput, MVBU_Staninput)) {
    stop("staninput must be an MVBU_Staninput object.", call. = FALSE)
  }

  x@staninput <- staninput
  x
}

S7::method(set_staninput, list(IdealAdaptorStanfitInput, S7::class_any)) <- function(x, staninput) {
  if (!S7::S7_inherits(staninput, MVBU_Staninput)) {
    stop("staninput must be an MVBU_Staninput object.", call. = FALSE)
  }

  x@staninput <- staninput
  x
}

S7::method(get_transform_information, IdealAdaptorStanfitInput) <- function(x) {
  x@transform_information
}
