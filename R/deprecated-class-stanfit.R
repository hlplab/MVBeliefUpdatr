#' Deprecated: is.ideal_adaptor_stanfit
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#' `is.ideal_adaptor_stanfit()` is deprecated; use
#' `S7::S7_inherits(x, IdealAdaptorStanfit)` instead.
#'
#' @param x Object to be checked.
#' @param verbose Currently being ignored.
#' @return A logical.
#' @keywords internal
#' @rdname deprecated-functions
#' @export
is.ideal_adaptor_stanfit <- function(x, verbose = FALSE) {
    lifecycle::deprecate_warn(
        "0.2.0",
        "is.ideal_adaptor_stanfit()",
        details = "Use S7::S7_inherits(x, IdealAdaptorStanfit) instead."
    )
    S7::S7_inherits(x, IdealAdaptorStanfit)
}

#' Deprecated: ideal_adaptor_stanfit
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#' `ideal_adaptor_stanfit()` is deprecated; use the S7 constructor
#' [MVBU_Stanfit()] or [fit_ideal_adaptor()] instead.
#'
#' @inheritParams MVBU_Stanfit
#' @return An object of class `MVBU_Stanfit`.
#' @keywords internal
#' @rdname deprecated-functions
#' @export
ideal_adaptor_stanfit <- function(
  data = data.frame(),
  staninput = NULL,
  stanvars = NULL,
  backend = "rstan",
  save_pars = NULL,
  stan_args = list(),
  stanfit = NULL,
  basis = NULL,
  transform_information = NULL,
  criteria = list(),
  file = NULL,
  version = NULL,
  metadata = list()
) {
  lifecycle::deprecate_warn(
    "0.1.0",
    "ideal_adaptor_stanfit()",
    details = "Use MVBU_Stanfit() or fit_ideal_adaptor() instead."
  )
  constructor <- .get_ideal_adaptor_stanfit_constructor(staninput)

  constructor(
    data = data,
    staninput = staninput,
    stanvars = stanvars,
    backend = backend,
    save_pars = save_pars,
    stan_args = stan_args,
    stanfit = stanfit,
    basis = basis,
    transform_information = transform_information,
    criteria = criteria,
    file = file,
    version = version,
    metadata = metadata
  )
}
