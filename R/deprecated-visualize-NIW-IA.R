# =============================================================================
# Deprecated NIW-IA Visualization Functions
# =============================================================================

#' @include S7-generics.R
#' @include S7-plot-methods.R
NULL

# -----------------------------------------------------------------------------
# deprecated
# -----------------------------------------------------------------------------

#' Deprecated: plot_expected_categorization_function_1D
#'
#' @description `r lifecycle::badge("deprecated")`
#' `plot_expected_categorization_function_1D()` was deprecated in MVBeliefUpdatr 0.1.0 and will be
#' removed in 0.2.0. Please use [plot_categorization_functions()] instead.
#'
#' @param x Model object.
#' @param ... Arguments passed to [plot_categorization_functions()].
#' @seealso [plot_categorization_functions()]
#' @keywords internal
#' @export
plot_expected_categorization_function_1D <- function(x, ...) {
  lifecycle::deprecate_warn(
    "0.1.0",
    "plot_expected_categorization_function_1D()",
    with = "plot_categorization_functions()"
  )
  plot_categorization_functions(x, ...)
}

#' Deprecated: plot_expected_categorization_function_2D
#'
#' @description `r lifecycle::badge("deprecated")`
#' `plot_expected_categorization_function_2D()` was deprecated in MVBeliefUpdatr 0.1.0 and will be
#' removed in 0.2.0. Please use [plot_categorization_functions()] instead.
#'
#' @param x Model object.
#' @param ... Arguments passed to [plot_categorization_functions()].
#' @seealso [plot_categorization_functions()]
#' @keywords internal
#' @export
plot_expected_categorization_function_2D <- function(x, ...) {
  lifecycle::deprecate_warn(
    "0.1.0",
    "plot_expected_categorization_function_2D()",
    with = "plot_categorization_functions()"
  )
  plot_categorization_functions(x, ...)
}

