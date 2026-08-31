# =============================================================================
# Deprecated Expected Categories Plotting Functions
# =============================================================================

#' @include S7-generics.R
#' @include S7-plot-methods.R
NULL

# -----------------------------------------------------------------------------
# deprecated
# -----------------------------------------------------------------------------

#' Deprecated: Plot functions
#'
#' @name deprecated-expected-categories-plots
#' @rdname deprecated-functions
#' @keywords internal
NULL

#' Deprecated: plot_expected_categories
#' @rdname deprecated-functions
#' @export
plot_expected_categories <- function(model, ...) {
  lifecycle::deprecate_warn(
    when = "0.2.0",
    what = "plot_expected_categories()",
    with = "plot_categories()"
  )
  plot_categories(model, ...)
}

#' Deprecated: plot_expected_categories_contour
#' @rdname deprecated-functions
#' @export
plot_expected_categories_contour <- function(model, ...) {
  lifecycle::deprecate_warn(
    when = "0.2.0",
    what = "plot_expected_categories_contour()",
    with = "plot_categories()"
  )
  plot_categories(model, aes = "contour", ...)
}

#' Deprecated: plot_expected_categories_density
#' @rdname deprecated-functions
#' @export
plot_expected_categories_density <- function(model, ...) {
  lifecycle::deprecate_warn(
    when = "0.2.0",
    what = "plot_expected_categories_density()",
    with = "plot_categories()"
  )
  plot_categories(model, aes = "fill", ...)
}

#' Deprecated: plot_expected_categories_contour2D
#' @rdname deprecated-functions
#' @export
plot_expected_categories_contour2D <- function(model, ...) {
  lifecycle::deprecate_warn(
    when = "0.2.0",
    what = "plot_expected_categories_contour2D()",
    with = "plot_categories()"
  )
  plot_categories(model, aes = "contour", ...)
}

#' Deprecated: plot_expected_categories_density1D
#' @rdname deprecated-functions
#' @export
plot_expected_categories_density1D <- function(model, ...) {
  lifecycle::deprecate_warn(
    when = "0.2.0",
    what = "plot_expected_categories_density1D()",
    with = "plot_categories()"
  )
  plot_categories(model, aes = "fill", ...)
}

#' Deprecated: plot_expected_categories_density2D
#' @rdname deprecated-functions
#' @export
plot_expected_categories_density2D <- function(model, ...) {
  lifecycle::deprecate_warn(
    when = "0.2.0",
    what = "plot_expected_categories_density2D()",
    with = "plot_categories()"
  )
  plot_categories(model, aes = "fill", ...)
}

#' Deprecated: plot_expected_categories.ideal_adaptor_stanfit
#' @rdname deprecated-functions
#' @export
plot_expected_categories.ideal_adaptor_stanfit <- function(
  model,
  type = "density",
  cues = get_cue_labels(model),
  ...
) {
  lifecycle::deprecate_warn(
    when = "0.2.0",
    what = "plot_expected_categories.ideal_adaptor_stanfit()",
    with = "plot_categories()"
  )
  if (type == "contour") {
    if (length(cues) == 1L) {
      warning(paste0(
        "Contour plots are only supported when at least 2 cues ",
        "are selected for plotting."
      ))
      return(plot_categories(model, cues = cues, aes = "fill", ...))
    }
    if (length(cues) > 2L) {
      warning("Contour plots for more than 2 are not yet supported.")
      return(invisible(NULL))
    }
    return(plot_categories(model, cues = cues, aes = "contour", ...))
  } else {
    if (length(cues) > 2L) {
      warning("Density plots for more than 2 are not yet supported.")
      return(invisible(NULL))
    }
    return(plot_categories(model, cues = cues, aes = "fill", ...))
  }
}
