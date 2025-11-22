#' Extract Diagnostic Quantities of \pkg{MVBeliefUpdatr} Models
#'
#' Extract quantities that can be used to diagnose sampling behavior
#' of the algorithms applied by \pkg{Stan} at the back-end of \pkg{MVBeliefUpdatr}.
#' These diagnostic functions are all copied and modified from \pkg{brms}.
#'
#' @name diagnostic-quantities
#' @aliases log_posterior nuts_params rhat neff_ratio
#'
#' @param x An \code{MVBeliefUpdatr} object.
#' @param pars An optional character vector of parameter names.
#'   For \code{nuts_params} these will be NUTS sampler parameter
#'   names rather than model parameters. If pars is omitted
#'   all parameters are included.
#' @param ... Arguments passed to individual methods.
#'
#' @return The exact form of the output depends on the method.
#'
#' @details For more details see
#'   \code{\link[bayesplot:bayesplot-extractors]{bayesplot-extractors}}.
#'
NULL

#' @rdname diagnostic-quantities
#' @importFrom bayesplot log_posterior
#' @export log_posterior
#' @export
log_posterior.ideal_adaptor_stanfit <- function(x, ...) {
  assert_contains_draws(x)
  bayesplot::log_posterior(x$stanfit, ...)
}

#' @rdname diagnostic-quantities
#' @importFrom bayesplot nuts_params
#' @export nuts_params
#' @export
nuts_params.ideal_adaptor_stanfit <- function(x, pars = NULL, ...) {
  assert_contains_draws(x)
  bayesplot::nuts_params(x$stanfit, pars = pars, ...)
}

#' @rdname diagnostic-quantities
#' @importFrom posterior rhat
#' @export rhat
#' @export
rhat.ideal_adaptor_stanfit <- function(x, pars = NULL, ...) {
  assert_contains_draws(x)
  # bayesplot uses outdated rhat code from rstan
  # bayesplot::rhat(x$stanfit, pars = pars, ...)
  draws <- as_draws_array(x, variable = pars, ...)
  tmp <- posterior::summarise_draws(draws, rhat = posterior::rhat)
  rhat <- tmp$rhat
  names(rhat) <- tmp$variable
  rhat
}

#' @rdname diagnostic-quantities
#' @importFrom bayesplot neff_ratio
#' @export neff_ratio
#' @export
neff_ratio.ideal_adaptor_stanfit <- function(x, pars = NULL, ...) {
  assert_contains_draws(x)
  # bayesplot uses outdated ess code from rstan
  # bayesplot::neff_ratio(x$stanfit, pars = pars, ...)
  draws <- as_draws_array(x, variable = pars, ...)
  tmp <- posterior::summarise_draws(
    draws, ess_bulk = posterior::ess_bulk, ess_tail = posterior::ess_tail
  )
  # min of ess_bulk and ess_tail mimics definition of posterior::rhat.default
  ess <- matrixStats::rowMins(cbind(tmp$ess_bulk, tmp$ess_tail))
  names(ess) <- tmp$variable
  ess / ndraws(draws)
}

#' Extract Control Parameters of the NUTS Sampler
#'
#' Extract control parameters of the NUTS sampler such as
#' \code{adapt_delta} or \code{max_treedepth}.
#'
#' @param x An \R object
#' @param pars Optional names of the control parameters to be returned.
#'  If \code{NULL} (the default) all control parameters are returned.
#'  See \code{\link[rstan]{stan}} for more details.
#' @param ... Currently ignored.
#'
#' @return A named \code{list} with control parameter values.
#'
#' @export
control_params <- function(x, ...) {
  UseMethod("control_params")
}

#' @rdname control_params
#' @export
control_params.ideal_adaptor_stanfit <- function(x, pars = NULL, ...) {
  assert_contains_draws(x)

  if (is_equal(x$backend, "cmdstanr")) {
    out <- attr(x$stanfit, "metadata")$metadata
  } else {
    out <- attr(x$stanfit@sim$samples[[1]], "args")$control
  }
  if (!is.null(pars)) {
    out <- out[pars]
  }
  out
}
