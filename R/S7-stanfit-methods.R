#' @include S7-core-classes.R
#' @include S7-generics.R
#' @include S7-transform-information.R
#' @include S7-staninput.R
#' @include S7-stanfit-input.R
#' @include S7-stanfit.R
NULL

S7::method(get_stanfit, S7::class_any) <- function(x) {
  stop("x must be an IdealAdaptorStanfit object.", call. = FALSE)
}

S7::method(get_stanfit, MVBU_Stanfit) <- function(x) {
  x@stanfit
}

S7::method(set_stanfit, list(S7::class_any, S7::class_any)) <- function(x, stanfit) {
  stop("x must be an IdealAdaptorStanfit object.", call. = FALSE)
}

S7::method(set_stanfit, list(MVBU_Stanfit, S7::class_any)) <- function(x, stanfit) {
  if (!is.null(stanfit)) {
    .assert_stanfit(stanfit)
    .assert_that(
      stanfit@model_name %in% names(MVBeliefUpdatr:::stanmodels),
      msg = paste0(
        "stanfit object was not created by one of the accepted stancodes:\n\t",
        paste(names(MVBeliefUpdatr:::stanmodels), collapse = "\n\t"),
        "\n(you can get the name of your model from your_stanfit@model_name)."
      )
    )
  }

  constructor <- get_ideal_adaptor_stanfit_constructor(x@staninput)

  constructor(
    data = x@data,
    staninput = x@staninput,
    stanvars = x@stanvars,
    backend = x@backend,
    save_pars = x@save_pars,
    stan_args = x@stan_args,
    stanfit = stanfit,
    basis = x@basis,
    transform_information = x@transform_information,
    criteria = x@criteria,
    file = x@file,
    version = x@version,
    labels = x@labels
  )
}

S7::method(get_staninput, S7::class_any) <- function(x) {
  stop("x must be an IdealAdaptorStanfit or IdealAdaptorStanfitInput object.", call. = FALSE)
}

S7::method(get_staninput, MVBU_Stanfit) <- function(x) {
  x@staninput
}

S7::method(get_transform_information, S7::class_any) <- function(x) {
  stop("x must be an IdealAdaptorStanfit or IdealAdaptorStanfitInput object.", call. = FALSE)
}

S7::method(get_transform_information, MVBU_Stanfit) <- function(x) {
  x@transform_information
}
