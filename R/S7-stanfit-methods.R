#' @include asserts.R
#' @include S7-core-classes.R
#' @include S7-generics.R
#' @include S7-transform-information.R
#' @include S7-staninput.R
#' @include S7-stanfit-input.R
#' @include S7-stanfit.R
NULL

get_ideal_adaptor_stanfit_constructor <- function(staninput = NULL) {
  if (!is.null(staninput)) {
    if (S7::S7_inherits(staninput, NIX_IdealAdaptorStaninput)) {
      NIX_IdealAdaptorStanfit
    } else if (S7::S7_inherits(staninput, MNIX_IdealAdaptorStaninput)) {
      MNIX_IdealAdaptorStanfit
    } else if (S7::S7_inherits(staninput, NIW_IdealAdaptorStaninput)) {
      NIW_IdealAdaptorStanfit
    } else {
      IdealAdaptorStanfit
    }
  } else {
    IdealAdaptorStanfit
  }
}

S7::method(get_stanfit, S7::class_any) <- function(x) {
  .stop("x must be an IdealAdaptorStanfit object.")
}

S7::method(get_stanfit, MVBU_Stanfit) <- function(x) {
  x@stanfit
}

S7::method(set_stanfit, list(S7::class_any, S7::class_any)) <- function(x, stanfit) {
  .stop("x must be an IdealAdaptorStanfit object.")
}

S7::method(set_stanfit, list(MVBU_Stanfit, S7::class_any)) <- function(x, stanfit) {
  # no assertions for stanfit here since the @<- assignment operator applied to S7 objects
  # will automatically call the validator for the class, which already check that the stanfit 
  # is valid.
  x@stanfit <- stanfit
  x
}

S7::method(get_staninput, S7::class_any) <- function(x) {
  .stop("x must be an IdealAdaptorStanfit or IdealAdaptorStanfitInput object.")
}

S7::method(get_staninput, MVBU_Stanfit) <- function(x) {
  x@staninput
}

S7::method(get_transform_information, S7::class_any) <- function(x) {
  .stop("x must be an IdealAdaptorStanfit or IdealAdaptorStanfitInput object.")
}

S7::method(get_transform_information, MVBU_Stanfit) <- function(x) {
  x@transform_information
}
