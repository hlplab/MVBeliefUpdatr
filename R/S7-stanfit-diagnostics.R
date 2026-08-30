#' @include asserts.R
#' @include S7-core-classes.R
#' @include S7-generics.R
#' @include S7-stanfit.R
#' @include S7-stanfit-methods.R
NULL

#' Extract Diagnostic Quantities of \pkg{MVBeliefUpdatr} Stanfit Models
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
#'   [bayesplot::bayesplot-extractors].
#'
#' @export log_posterior
#' @export nuts_params
#' @export rhat
#' @export neff_ratio
NULL

#' @rdname diagnostic-quantities
#' @export
S7::method(log_posterior, MVBU_Stanfit) <- function(x, ...) {
  .assert_contains_draws(x)
  bayesplot::log_posterior(get_stanfit(x), ...)
}

S7::method(log_posterior, S7::class_any) <- function(x, ...) {
  bayesplot::log_posterior(x, ...)
}

#' @rdname diagnostic-quantities
#' @export
S7::method(nuts_params, MVBU_Stanfit) <- function(x, pars = NULL, ...) {
  .assert_contains_draws(x)
  bayesplot::nuts_params(get_stanfit(x), pars = pars, ...)
}

S7::method(nuts_params, S7::class_any) <- function(x, pars = NULL, ...) {
  bayesplot::nuts_params(x, pars = pars, ...)
}

#' @rdname diagnostic-quantities
#' @export
S7::method(rhat, MVBU_Stanfit) <- function(x, pars = NULL, ...) {
  .assert_contains_draws(x)
  draws <- posterior::as_draws_array(get_stanfit(x), variable = pars, ...)
  tmp <- posterior::summarise_draws(draws, rhat = posterior::rhat)
  rhat_vals <- tmp$rhat
  names(rhat_vals) <- tmp$variable
  rhat_vals
}

S7::method(rhat, S7::class_any) <- function(x, ...) {
  posterior::rhat(x, ...)
}

#' @rdname diagnostic-quantities
#' @export
S7::method(neff_ratio, MVBU_Stanfit) <- function(x, pars = NULL, ...) {
  .assert_contains_draws(x)
  draws <- posterior::as_draws_array(get_stanfit(x), variable = pars, ...)
  tmp <- posterior::summarise_draws(
    draws,
    ess_bulk = posterior::ess_bulk, ess_tail = posterior::ess_tail
  )
  # min of ess_bulk and ess_tail mimics definition of posterior::rhat.default
  ess <- matrixStats::rowMins(cbind(tmp$ess_bulk, tmp$ess_tail))
  names(ess) <- tmp$variable
  ess / posterior::ndraws(draws)
}

S7::method(neff_ratio, S7::class_any) <- function(x, pars = NULL, ...) {
  bayesplot::neff_ratio(x, pars = pars, ...)
}

#' @rdname control_params
#' @export
S7::method(control_params, MVBU_Stanfit) <- function(x, pars = NULL, ...) {
  .assert_contains_draws(x)
  backend <- x@backend
  sf <- get_stanfit(x)
  if (.is_equal(backend, "cmdstanr")) {
    out <- attr(sf, "metadata")$metadata
  } else {
    out <- attr(sf@sim$samples[[1]], "args")$control
  }
  if (!is.null(pars)) {
    out <- out[pars]
  }
  out
}

S7::method(control_params, S7::class_any) <- function(x, pars = NULL, ...) {
  if (inherits(x, "stanfit")) {
    out <- attr(x@sim$samples[[1]], "args")$control
    if (!is.null(pars)) out <- out[pars]
    out
  } else {
    .stop("control_params is not implemented for objects of class ", class(x)[1])
  }
}

# Kept as_draws*.MVBU_Stanfit S3 methods for the external posterior package so that posterior::as_draws_df(fit)
# and related functions dispatch as expected when called directly on MVBU_Stanfit.

#' @rdname diagnostic-quantities
#' @importFrom posterior as_draws as_draws_df as_draws_array as_draws_matrix as_draws_list as_draws_rvars
#' @export
as_draws.MVBU_Stanfit <- function(x, ...) {
  posterior::as_draws(get_stanfit(x), ...)
}

#' @rdname diagnostic-quantities
#' @export
as_draws_df.MVBU_Stanfit <- function(x, ...) {
  posterior::as_draws_df(get_stanfit(x), ...)
}

#' @rdname diagnostic-quantities
#' @export
as_draws_array.MVBU_Stanfit <- function(x, ...) {
  posterior::as_draws_array(get_stanfit(x), ...)
}

#' @rdname diagnostic-quantities
#' @export
as_draws_matrix.MVBU_Stanfit <- function(x, ...) {
  posterior::as_draws_matrix(get_stanfit(x), ...)
}

#' @rdname diagnostic-quantities
#' @export
as_draws_list.MVBU_Stanfit <- function(x, ...) {
  posterior::as_draws_list(get_stanfit(x), ...)
}

#' @rdname diagnostic-quantities
#' @export
as_draws_rvars.MVBU_Stanfit <- function(x, ...) {
  posterior::as_draws_rvars(get_stanfit(x), ...)
}

