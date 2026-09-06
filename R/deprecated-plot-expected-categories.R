# =============================================================================
# Deprecated Expected Categories Plotting Functions
# =============================================================================

#' @include S7-generics.R
#' @include S7-plot-methods.R
NULL

# -----------------------------------------------------------------------------
# deprecated
# -----------------------------------------------------------------------------

#' Deprecated: plot_expected_categories
#'
#' @description `r lifecycle::badge("deprecated")`
#' `plot_expected_categories()` was deprecated in MVBeliefUpdatr 0.1.0 and will be
#' removed in 0.2.0. Please use [plot_categories()] instead.
#'
#' @param model Model object.
#' @param ... Arguments passed to [plot_categories()].
#' @seealso [plot_categories()]
#' @keywords internal
#' @export
plot_expected_categories <- function(model, ...) {
  lifecycle::deprecate_warn(
    "0.1.0",
    "plot_expected_categories()",
    with = "plot_categories()"
  )
  plot_categories(model, ...)
}

#' Deprecated: plot_expected_categories_contour
#'
#' @description `r lifecycle::badge("deprecated")`
#' `plot_expected_categories_contour()` was deprecated in MVBeliefUpdatr 0.1.0 and will be
#' removed in 0.2.0. Please use [plot_categories()] instead.
#'
#' @param model Model object.
#' @param ... Arguments passed to [plot_categories()].
#' @seealso [plot_categories()]
#' @keywords internal
#' @export
plot_expected_categories_contour <- function(model, ...) {
  lifecycle::deprecate_warn(
    "0.1.0",
    "plot_expected_categories_contour()",
    with = "plot_categories()"
  )
  plot_categories(model, aes = "contour", ...)
}

#' Deprecated: plot_expected_categories_density
#'
#' @description `r lifecycle::badge("deprecated")`
#' `plot_expected_categories_density()` was deprecated in MVBeliefUpdatr 0.1.0 and will be
#' removed in 0.2.0. Please use [plot_categories()] instead.
#'
#' @param model Model object.
#' @param ... Arguments passed to [plot_categories()].
#' @seealso [plot_categories()]
#' @keywords internal
#' @export
plot_expected_categories_density <- function(model, ...) {
  lifecycle::deprecate_warn(
    "0.1.0",
    "plot_expected_categories_density()",
    with = "plot_categories()"
  )
  plot_categories(model, aes = "fill", ...)
}

#' Deprecated: plot_expected_categories_contour2D
#'
#' @description `r lifecycle::badge("deprecated")`
#' `plot_expected_categories_contour2D()` was deprecated in MVBeliefUpdatr 0.1.0 and will be
#' removed in 0.2.0. Please use [plot_categories()] instead.
#'
#' @param model Model object.
#' @param ... Arguments passed to [plot_categories()].
#' @seealso [plot_categories()]
#' @keywords internal
#' @export
plot_expected_categories_contour2D <- function(model, ...) {
  lifecycle::deprecate_warn(
    "0.1.0",
    "plot_expected_categories_contour2D()",
    with = "plot_categories()"
  )
  plot_categories(model, aes = "contour", ...)
}

#' Deprecated: plot_expected_categories_density1D
#'
#' @description `r lifecycle::badge("deprecated")`
#' `plot_expected_categories_density1D()` was deprecated in MVBeliefUpdatr 0.1.0 and will be
#' removed in 0.2.0. Please use [plot_categories()] instead.
#'
#' @param model Model object.
#' @param ... Arguments passed to [plot_categories()].
#' @seealso [plot_categories()]
#' @keywords internal
#' @export
plot_expected_categories_density1D <- function(model, ...) {
  lifecycle::deprecate_warn(
    "0.1.0",
    "plot_expected_categories_density1D()",
    with = "plot_categories()"
  )
  plot_categories(model, aes = "fill", ...)
}

#' Deprecated: plot_expected_categories_density2D
#'
#' @description `r lifecycle::badge("deprecated")`
#' `plot_expected_categories_density2D()` was deprecated in MVBeliefUpdatr 0.1.0 and will be
#' removed in 0.2.0. Please use [plot_categories()] instead.
#'
#' @param model Model object.
#' @param ... Arguments passed to [plot_categories()].
#' @seealso [plot_categories()]
#' @keywords internal
#' @export
plot_expected_categories_density2D <- function(model, ...) {
  lifecycle::deprecate_warn(
    "0.1.0",
    "plot_expected_categories_density2D()",
    with = "plot_categories()"
  )
  plot_categories(model, aes = "fill", ...)
}

#' Deprecated: plot_expected_categories.ideal_adaptor_stanfit
#'
#' @description `r lifecycle::badge("deprecated")`
#' `plot_expected_categories.ideal_adaptor_stanfit()` was deprecated in MVBeliefUpdatr 0.1.0 and will be
#' removed in 0.2.0. Please use [plot_categories()] instead.
#'
#' @param model Model object.
#' @param type Plot type ("density" or "contour").
#' @param cues Cue names to plot.
#' @param ... Arguments passed to [plot_categories()].
#' @seealso [plot_categories()]
#' @keywords internal
#' @export
plot_expected_categories.ideal_adaptor_stanfit <- function(
  model,
  type = "density",
  cues = get_cue_labels(model),
  ...
) {
  lifecycle::deprecate_warn(
    "0.1.0",
    "plot_expected_categories.ideal_adaptor_stanfit()",
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

