#' @include S7-core-classes.R
#' @include S7-generics.R
#' @include S7-stanfit.R
#' @include S7-stanfit-methods.R
#' @importFrom tidybayes recover_types
#' @importFrom lifecycle deprecate_warn
#' @importFrom tibble tibble
#' @importFrom dplyr group_by summarise mutate
#' @importFrom purrr reduce
#' @importFrom rlang sym syms .data :=
NULL

# deprecated ------------------------------------------------------------------


#' Get the name of the stanmodel from an ideal adaptor stanfit (deprecated)
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#' `get_stanmodelname()` is deprecated; use \code{\link{get_model_type}} instead.
#'
#' @param x A model object.
#' @param ... Additional arguments.
#' @return A character string.
#' @export
get_stanmodelname <- function(x, ...) {
  lifecycle::deprecate_warn(
    "0.0.9",
    "get_stanmodelname()",
    "get_model_type()"
  )
  get_model_type(x)
}

#' Get or restore the original group or category levels (deprecated)
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#' `get_staninput_variable_levels()` is deprecated in favor of
#' \code{\link{get_category_labels}}, \code{\link{get_group_labels}}, and
#' \code{\link{get_cue_labels}}.
#'
#' @param x A model object.
#' @param variable Either "category", "group", or "cue".
#' @param indices Optional indices.
#' @export
get_staninput_variable_levels <- function(
  x,
  variable = c("category", "group", "cue"),
  indices = NULL
) {
  lifecycle::deprecate_warn(
    "0.0.9",
    "get_staninput_variable_levels()",
    "get_labels()"
  )
  variable <- match.arg(variable)
  switch(variable,
    "category" = get_category_labels(x, indices = indices),
    "group" = get_group_labels(x, indices = indices),
    "cue" = get_cue_labels(x, indices = indices)
  )
}

#' Get category levels from model (deprecated)
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#' `get_category_levels()` is deprecated in favor of
#' \code{\link{get_category_labels}}.
#'
#' @param x A model object.
#' @param indices Optional indices.
#' @param ... Additional arguments.
#' @export
get_category_levels <- function(x, indices = NULL, ...) {
  lifecycle::deprecate_warn(
    "0.0.9",
    "get_category_levels()",
    "get_category_labels()"
  )
  get_category_labels(x, indices = indices, ...)
}

#' Get group levels from model (deprecated)
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#' `get_group_levels()` is deprecated in favor of
#' \code{\link{get_group_labels}}.
#'
#' @param x A model object.
#' @param indices Optional indices.
#' @param include_prior Whether to include `"prior"`.
#' @param ... Additional arguments.
#' @export
get_group_levels <- function(x, indices = NULL, include_prior = FALSE, ...) {
  lifecycle::deprecate_warn(
    "0.0.9",
    "get_group_levels()",
    "get_group_labels()"
  )
  get_group_labels(x, indices = indices, include_prior = include_prior, ...)
}

#' Get cue levels from model (deprecated)
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#' `get_cue_levels()` is deprecated in favor of
#' \code{\link{get_cue_labels}}.
#'
#' @param x A model object.
#' @param indices Optional indices.
#' @param ... Additional arguments.
#' @export
get_cue_levels <- function(x, indices = NULL, ...) {
  lifecycle::deprecate_warn(
    "0.0.9",
    "get_cue_levels()",
    "get_cue_labels()"
  )
  get_cue_labels(x, indices = indices, ...)
}

#' Get parameter names from stanfit (deprecated)
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#' `get_params()` is deprecated in favor of
#' \code{\link{get_parameter_names}}.
#'
#' @param fit A stanfit or MVBU_Stanfit object.
#' @param original_pars Whether to return original parameter names.
#' @export
get_params <- function(fit, original_pars = FALSE) {
  lifecycle::deprecate_warn(
    "0.0.9",
    "get_params()",
    "get_parameter_names()"
  )
  get_parameter_names(fit, original_pars = original_pars)
}
