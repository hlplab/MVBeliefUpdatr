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


# deprecated ------------------------------------------------------------------


#' Deprecated: get_stanmodelname
#'
#' @description `r lifecycle::badge("deprecated")`
#' `get_stanmodelname()` was deprecated in MVBeliefUpdatr 0.1.0 and will be
#' removed in 0.2.0. Please use [get_model_type()] instead.
#'
#' @param x A model object.
#' @param ... Additional arguments.
#' @return A character string.
#' @seealso [get_model_type()]
#' @keywords internal
#' @export
get_stanmodelname <- function(x, ...) {
  lifecycle::deprecate_warn(
    "0.1.0",
    "get_stanmodelname()",
    with = "get_model_type()"
  )
  get_model_type(x)
}

#' Deprecated: get_staninput_variable_levels
#'
#' @description `r lifecycle::badge("deprecated")`
#' `get_staninput_variable_levels()` was deprecated in MVBeliefUpdatr 0.1.0 and will be
#' removed in 0.2.0. Please use [get_category_labels()], [get_group_labels()], or [get_cue_labels()] instead.
#'
#' @param x A model object.
#' @param variable Either "category", "group", or "cue".
#' @param indices Optional indices.
#' @seealso [get_category_labels()], [get_group_labels()], [get_cue_labels()]
#' @keywords internal
#' @export
get_staninput_variable_levels <- function(
  x,
  variable = c("category", "group", "cue"),
  indices = NULL
) {
  lifecycle::deprecate_warn(
    "0.1.0",
    "get_staninput_variable_levels()",
    details = "Use get_category_labels(), get_group_labels(), or get_cue_labels() instead."
  )
  variable <- match.arg(variable)
  switch(variable,
    "category" = get_category_labels(x, indices = indices),
    "group" = get_group_labels(x, indices = indices),
    "cue" = get_cue_labels(x, indices = indices)
  )
}

#' Deprecated: get_category_levels
#'
#' @description `r lifecycle::badge("deprecated")`
#' `get_category_levels()` was deprecated in MVBeliefUpdatr 0.1.0 and will be
#' removed in 0.2.0. Please use [get_category_labels()] instead.
#'
#' @param x A model object.
#' @param indices Optional indices.
#' @param ... Additional arguments.
#' @seealso [get_category_labels()]
#' @keywords internal
#' @export
get_category_levels <- function(x, indices = NULL, ...) {
  lifecycle::deprecate_warn(
    "0.1.0",
    "get_category_levels()",
    with = "get_category_labels()"
  )
  get_category_labels(x, indices = indices, ...)
}

#' Deprecated: get_group_levels
#'
#' @description `r lifecycle::badge("deprecated")`
#' `get_group_levels()` was deprecated in MVBeliefUpdatr 0.1.0 and will be
#' removed in 0.2.0. Please use [get_group_labels()] instead.
#'
#' @param x A model object.
#' @param indices Optional indices.
#' @param include_prior Whether to include `"prior"`.
#' @param ... Additional arguments.
#' @seealso [get_group_labels()]
#' @keywords internal
#' @export
get_group_levels <- function(x, indices = NULL, include_prior = FALSE, ...) {
  lifecycle::deprecate_warn(
    "0.1.0",
    "get_group_levels()",
    with = "get_group_labels()"
  )
  get_group_labels(x, indices = indices, include_prior = include_prior, ...)
}

#' Deprecated: get_cue_levels
#'
#' @description `r lifecycle::badge("deprecated")`
#' `get_cue_levels()` was deprecated in MVBeliefUpdatr 0.1.0 and will be
#' removed in 0.2.0. Please use [get_cue_labels()] instead.
#'
#' @param x A model object.
#' @param indices Optional indices.
#' @param ... Additional arguments.
#' @seealso [get_cue_labels()]
#' @keywords internal
#' @export
get_cue_levels <- function(x, indices = NULL, ...) {
  lifecycle::deprecate_warn(
    "0.1.0",
    "get_cue_levels()",
    with = "get_cue_labels()"
  )
  get_cue_labels(x, indices = indices, ...)
}

#' Deprecated: get_params
#'
#' @description `r lifecycle::badge("deprecated")`
#' `get_params()` was deprecated in MVBeliefUpdatr 0.1.0 and will be
#' removed in 0.2.0. Please use [get_parameter_names()] instead.
#'
#' @param fit A stanfit or MVBU_Stanfit object.
#' @param original_pars Whether to return original parameter names.
#' @seealso [get_parameter_names()]
#' @keywords internal
#' @export
get_params <- function(fit, original_pars = FALSE) {
  lifecycle::deprecate_warn(
    "0.1.0",
    "get_params()",
    with = "get_parameter_names()"
  )
  get_parameter_names(fit, original_pars = original_pars)
}

