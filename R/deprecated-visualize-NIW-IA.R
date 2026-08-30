# =============================================================================
# Deprecated NIW-IA Visualization Functions
# =============================================================================

#' @include S7-generics.R
#' @include S7-plot-methods.R
NULL

# -----------------------------------------------------------------------------
# deprecated
# -----------------------------------------------------------------------------

#' Deprecated NIW-IA Visualization Functions
#'
#' @name deprecated-niw-ia-plots
#' @rdname deprecated-functions
#' @keywords internal
NULL

#' @rdname deprecated-functions
#' @export
plot_expected_categorization_function_1D <- function(x, ...) {
  lifecycle::deprecate_warn(
    when = "0.2.0",
    what = "plot_expected_categorization_function_1D()",
    with = "plot_categorization_function()"
  )
  plot_categorization_function(x, ...)
}

#' @rdname deprecated-functions
#' @export
plot_expected_categorization_function_2D <- function(x, ...) {
  lifecycle::deprecate_warn(
    when = "0.2.0",
    what = "plot_expected_categorization_function_2D()",
    with = "plot_categorization_function()"
  )
  plot_categorization_function(x, ...)
}
