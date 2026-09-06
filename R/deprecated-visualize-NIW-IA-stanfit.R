# =============================================================================
# Deprecated Stanfit Plotting Functions
# =============================================================================

#' @include S7-generics.R
#' @include S7-plot-methods.R
NULL

# -----------------------------------------------------------------------------
# deprecated
# -----------------------------------------------------------------------------

#' Deprecated: plot_parameters.ideal_adaptor_stanfit
#'
#' @description `r lifecycle::badge("deprecated")`
#' `plot_parameters.ideal_adaptor_stanfit()` was deprecated in MVBeliefUpdatr 0.1.0 and will be
#' removed in 0.2.0. Please use [plot_parameters()] instead.
#'
#' @param model Model object.
#' @param ... Arguments passed to [plot_parameters()].
#' @seealso [plot_parameters()]
#' @keywords internal
#' @export
plot_parameters.ideal_adaptor_stanfit <- function(model, ...) {
  lifecycle::deprecate_warn(
    "0.1.0",
    "plot_parameters.ideal_adaptor_stanfit()",
    with = "plot_parameters()"
  )
  plot_parameters(model, ...)
}

#' Deprecated: plot_parameter_correlations.ideal_adaptor_stanfit
#'
#' @description `r lifecycle::badge("deprecated")`
#' `plot_parameter_correlations.ideal_adaptor_stanfit()` was deprecated in MVBeliefUpdatr 0.1.0 and will be
#' removed in 0.2.0. Please use [plot_parameter_correlations()] instead.
#'
#' @param model Model object.
#' @param ... Arguments passed to [plot_parameter_correlations()].
#' @seealso [plot_parameter_correlations()]
#' @keywords internal
#' @export
plot_parameter_correlations.ideal_adaptor_stanfit <- function(model, ...) {
  lifecycle::deprecate_warn(
    "0.1.0",
    "plot_parameter_correlations.ideal_adaptor_stanfit()",
    with = "plot_parameter_correlations()"
  )
  plot_parameter_correlations(model, ...)
}

#' Deprecated: plot_expected_categorization_function_from_stanfit
#'
#' @description `r lifecycle::badge("deprecated")`
#' `plot_expected_categorization_function_from_stanfit()` was deprecated in MVBeliefUpdatr 0.1.0 and will be
#' removed in 0.2.0. Please use [plot_categorization_function()] instead.
#'
#' @param model Model object.
#' @param ... Arguments passed to [plot_categorization_function()].
#' @seealso [plot_categorization_function()]
#' @keywords internal
#' @export
plot_expected_categorization_function_from_stanfit <- function(
  model,
  ...
) {
  lifecycle::deprecate_warn(
    "0.1.0",
    "plot_expected_categorization_function_from_stanfit()",
    with = "plot_categorization_function()"
  )
  plot_categorization_function(model, ...)
}

