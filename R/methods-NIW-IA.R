#' Summarize NIW ideal adaptor stanfit
#'
#' \code{summary} method for \code{\link{NIW_ideal_adaptor_stanfit}} objects. Specifies reasonable defaults for the parameters to be summarized for the stanfit object.
#'
#' @param x An \code{\link{NIW_ideal_adaptor_stanfit}} object.
#'
#' @method summary NIW_ideal_adaptor_stanfit
#'
#' @export
#' @export summary
#' @importFrom rstan summary
summary.NIW_ideal_adaptor_stanfit <- function(x, pars = NULL, ...) {
  if (is.null(pars)) {
    pars <- names(x)
    pars <- grep("^((kappa|nu|m|S)_|lapse_rate)", pars, value = T)
    pars <- grep("^m_0_(tau|L_omega)", pars, value = T, invert = T)
    pars <- grep("^(m|S)_0_param", pars, value = T, invert = T)
  }

  rstan::summary(x, pars = pars, ...)$summary
}

#' loo NIW ideal adaptor stanfit
#'
#' \code{loo} method for \code{\link{NIW_ideal_adaptor_stanfit}} objects.
#'
#' @param x An \code{\link{NIW_ideal_adaptor_stanfit}} object.
#'
#' @method loo NIW_ideal_adaptor_stanfit
#'
#' @export loo
#' @importFrom loo loo extract_log_lik relative_eff loo.array
loo.NIW_ideal_adaptor_stanfit <- function(
    x,
    pars = "log_lik",
    ...,
    save_psis = FALSE,
    cores = getOption("mc.cores", 1)
) {
  stopifnot(length(pars) == 1L)
  LLarray <- loo::extract_log_lik(stanfit = x,
                                  parameter_name = pars,
                                  merge_chains = FALSE)
  r_eff <- loo::relative_eff(x = exp(LLarray), cores = cores)
  loo::loo.array(LLarray,
                 r_eff = r_eff,
                 cores = cores,
                 save_psis = save_psis)
}
