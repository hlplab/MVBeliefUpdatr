#' Deprecated: make_staninput
#'
#' @description `r lifecycle::badge("deprecated")`
#' `make_staninput()` was deprecated in MVBeliefUpdatr 0.1.0 and will be
#' removed in 0.2.0. Please use [new_ideal_adaptor_stanfit_input()] instead.
#'
#' @inheritParams new_ideal_adaptor_stanfit_input
#' @param tau_scale `r lifecycle::badge("deprecated")` Deprecated. Use
#'   `priors = prior(..., class = "tau_0")` in \code{\link{fit_ideal_adaptor}}
#'   instead.
#' @param L_omega_eta `r lifecycle::badge("deprecated")` Deprecated. Use
#'   `priors = prior(..., class = "cor_0")` in \code{\link{fit_ideal_adaptor}}
#'   instead.
#' @return An object of class \code{\link{IdealAdaptorStanfitInput}}.
#' @seealso [new_ideal_adaptor_stanfit_input()]
#' @keywords internal
#' @export
make_staninput <- function(
    exposure, test,
    cues, category = "category", response = "response",
    group = "group", group.unique = NULL,
    lapse_rate = NULL, mu_0 = NULL, Sigma_0 = NULL,
    tau_scale = NULL, L_omega_eta = NULL,
    control = control_staninput(),
    stanmodel = "NIW_ideal_adaptor",
    verbose = FALSE
) {
  lifecycle::deprecate_warn(
    "0.1.0",
    "make_staninput()",
    with = "new_ideal_adaptor_stanfit_input()"
  )
  if (!is.null(tau_scale)) {
    lifecycle::deprecate_warn(
      "0.9.0",
      "make_staninput(tau_scale)",
      details = paste0(
        "Use priors = prior(..., class = 'tau_0') in fit_ideal_adaptor() ",
        "instead."
      )
    )
  }
  if (!is.null(L_omega_eta)) {
    lifecycle::deprecate_warn(
      "0.9.0",
      "make_staninput(L_omega_eta)",
      details = paste0(
        "Use priors = prior(..., class = 'cor_0') in fit_ideal_adaptor() ",
        "instead."
      )
    )
  }
  if (!is.null(lapse_rate)) {
    lifecycle::deprecate_warn(
      "0.9.0",
      "make_staninput(lapse_rate)",
      details = paste0(
        "Use fixed_parameters = fix_parameter(lapse_rate = ...) in ",
        "fit_ideal_adaptor() instead."
      )
    )
  }
  if (!is.null(mu_0)) {
    lifecycle::deprecate_warn(
      "0.9.0",
      "make_staninput(mu_0)",
      details = paste0(
        "Use fixed_parameters = fix_parameter(m_0 = ...) in ",
        "fit_ideal_adaptor() instead."
      )
    )
  }
  if (!is.null(Sigma_0)) {
    lifecycle::deprecate_warn(
      "0.9.0",
      "make_staninput(Sigma_0)",
      details = paste0(
        "Use fixed_parameters = fix_parameter(expected_Sigma_0 = ...) in ",
        "fit_ideal_adaptor() instead."
      )
    )
  }
  new_ideal_adaptor_stanfit_input(
    exposure = exposure,
    test = test,
    cues = cues,
    category = category,
    response = response,
    group = group,
    group_unique = group.unique,
    control = control,
    stanmodel = stanmodel
  )
}

#' Deprecated: make_ideal_adaptor_stanfit_input
#'
#' @description `r lifecycle::badge("deprecated")`
#' `make_ideal_adaptor_stanfit_input()` was deprecated in MVBeliefUpdatr 0.1.0
#' and will be removed in 0.2.0. Please use [new_ideal_adaptor_stanfit_input()]
#' instead.
#'
#' @inheritParams new_ideal_adaptor_stanfit_input
#' @param tau_scale `r lifecycle::badge("deprecated")` Deprecated. Use
#'   `priors = prior(..., class = "tau_0")` in \code{\link{fit_ideal_adaptor}}
#'   instead.
#' @param L_omega_eta `r lifecycle::badge("deprecated")` Deprecated. Use
#'   `priors = prior(..., class = "cor_0")` in \code{\link{fit_ideal_adaptor}}
#'   instead.
#' @return An object of class \code{\link{IdealAdaptorStanfitInput}}.
#' @seealso [new_ideal_adaptor_stanfit_input()]
#' @keywords internal
#' @rdname make_staninput
#' @export
make_ideal_adaptor_stanfit_input <- function(
    exposure, test,
    cues, category = "category", response = "response",
    group = "group", group.unique = NULL,
    lapse_rate = NULL, mu_0 = NULL, Sigma_0 = NULL,
    tau_scale = NULL, L_omega_eta = NULL,
    control = control_staninput(),
    stanmodel = "NIW_ideal_adaptor",
    verbose = FALSE
) {
  lifecycle::deprecate_warn(
    "0.1.0",
    "make_ideal_adaptor_stanfit_input()",
    with = "new_ideal_adaptor_stanfit_input()"
  )
  if (!is.null(tau_scale)) {
    lifecycle::deprecate_warn(
      "0.9.0",
      "make_ideal_adaptor_stanfit_input(tau_scale)",
      details = paste0(
        "Use priors = prior(..., class = 'tau_0') in fit_ideal_adaptor() ",
        "instead."
      )
    )
  }
  if (!is.null(L_omega_eta)) {
    lifecycle::deprecate_warn(
      "0.9.0",
      "make_ideal_adaptor_stanfit_input(L_omega_eta)",
      details = paste0(
        "Use priors = prior(..., class = 'cor_0') in fit_ideal_adaptor() ",
        "instead."
      )
    )
  }
  if (!is.null(lapse_rate)) {
    lifecycle::deprecate_warn(
      "0.9.0",
      "make_ideal_adaptor_stanfit_input(lapse_rate)",
      details = paste0(
        "Use fixed_parameters = fix_parameter(lapse_rate = ...) in ",
        "fit_ideal_adaptor() instead."
      )
    )
  }
  if (!is.null(mu_0)) {
    lifecycle::deprecate_warn(
      "0.9.0",
      "make_ideal_adaptor_stanfit_input(mu_0)",
      details = paste0(
        "Use fixed_parameters = fix_parameter(m_0 = ...) in ",
        "fit_ideal_adaptor() instead."
      )
    )
  }
  if (!is.null(Sigma_0)) {
    lifecycle::deprecate_warn(
      "0.9.0",
      "make_ideal_adaptor_stanfit_input(Sigma_0)",
      details = paste0(
        "Use fixed_parameters = fix_parameter(expected_Sigma_0 = ...) in ",
        "fit_ideal_adaptor() instead."
      )
    )
  }
  new_ideal_adaptor_stanfit_input(
    exposure = exposure,
    test = test,
    cues = cues,
    category = category,
    response = response,
    group = group,
    group_unique = group.unique,
    control = control,
    stanmodel = stanmodel
  )
}
