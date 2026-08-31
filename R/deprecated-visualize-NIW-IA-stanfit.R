# =============================================================================
# Deprecated Stanfit Plotting Functions
# =============================================================================

#' @include S7-generics.R
#' @include S7-plot-methods.R
NULL

# -----------------------------------------------------------------------------
# deprecated
# -----------------------------------------------------------------------------

#' Deprecated Stanfit Plotting Functions
#'
#' @name deprecated-stanfit-plots
#' @rdname deprecated-functions
#' @keywords internal
NULL

#' Deprecated: plot_parameters.ideal_adaptor_stanfit
#' @rdname deprecated-functions
#' @export
plot_parameters.ideal_adaptor_stanfit <- function(model, ...) {
  lifecycle::deprecate_warn(
    when = "0.2.0",
    what = "plot_parameters.ideal_adaptor_stanfit()",
    with = "plot_parameters()"
  )
  plot_parameters(model, ...)
}

#' Deprecated: plot_parameter_correlations.ideal_adaptor_stanfit
#' @rdname deprecated-functions
#' @export
plot_parameter_correlations.ideal_adaptor_stanfit <- function(model, ...) {
  lifecycle::deprecate_warn(
    when = "0.2.0",
    what = "plot_parameter_correlations.ideal_adaptor_stanfit()",
    with = "plot_parameter_correlations()"
  )
  plot_parameter_correlations(model, ...)
}


#' Deprecated: plot_expected_categorization_function_from_stanfit
#' @rdname deprecated-functions
#' @export
plot_expected_categorization_function_from_stanfit <- function(
  model,
  ...
) {
  lifecycle::deprecate_warn(
    when = "0.2.0",
    what = "plot_expected_categorization_function_from_stanfit()",
    with = "plot_categorization_function()"
  )
  plot_categorization_function(model, ...)
}
