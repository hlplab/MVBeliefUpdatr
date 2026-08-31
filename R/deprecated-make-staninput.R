#' Deprecated: make_staninput
#'
#' @description Deprecated. Use \code{\link{new_ideal_adaptor_stanfit_input}} instead.
#' @inheritParams new_ideal_adaptor_stanfit_input
#' @return A list with components \code{staninput}, \code{data}, and \code{transform_information}.
#' @seealso \code{\link{new_ideal_adaptor_stanfit_input}}
#' @keywords internal
#' @export
make_staninput <- function(
    exposure, test,
    cues, category = "category", response = "response",
    group = "group", group.unique = NULL,
    lapse_rate = NULL, mu_0 = NULL, Sigma_0 = NULL,
    control = control_staninput(),
    stanmodel = "NIW_ideal_adaptor",
    verbose = FALSE
) {
  lifecycle::deprecate_warn(
    when = "0.0.3",
    what = "make_staninput()",
    with = "new_ideal_adaptor_stanfit_input()"
  )
  fixed_parameters <- list(
    lapse_rate = lapse_rate,
    mu_0 = mu_0,
    Sigma_0 = Sigma_0
  )
  new_ideal_adaptor_stanfit_input(
    exposure = exposure,
    test = test,
    cues = cues,
    category = category,
    response = response,
    group = group,
    group.unique = group.unique,
    fixed_parameters = fixed_parameters,
    control = control,
    stanmodel = stanmodel,
    verbose = verbose
  )
}

#' Deprecated: make_ideal_adaptor_stanfit_input
#'
#' @description Deprecated. Use \code{\link{new_ideal_adaptor_stanfit_input}} instead.
#' @inheritParams new_ideal_adaptor_stanfit_input
#' @return A list with components \code{staninput}, \code{data}, and \code{transform_information}.
#' @seealso \code{\link{new_ideal_adaptor_stanfit_input}}
#' @keywords internal
#' @rdname make_staninput
#' @export
make_ideal_adaptor_stanfit_input <- function(
    exposure, test,
    cues, category = "category", response = "response",
    group = "group", group.unique = NULL,
    lapse_rate = NULL, mu_0 = NULL, Sigma_0 = NULL,
    control = control_staninput(),
    stanmodel = "NIW_ideal_adaptor",
    verbose = FALSE
) {
  lifecycle::deprecate_warn(
    when = "0.0.3",
    what = "make_ideal_adaptor_stanfit_input()",
    with = "new_ideal_adaptor_stanfit_input()"
  )
  fixed_parameters <- list(
    lapse_rate = lapse_rate,
    mu_0 = mu_0,
    Sigma_0 = Sigma_0
  )
  new_ideal_adaptor_stanfit_input(
    exposure = exposure,
    test = test,
    cues = cues,
    category = category,
    response = response,
    group = group,
    group.unique = group.unique,
    fixed_parameters = fixed_parameters,
    control = control,
    stanmodel = stanmodel,
    verbose = verbose
  )
}
