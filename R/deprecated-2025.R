#' Deprecated: get_transform_information_from_stanfit
#'
#' @description `r lifecycle::badge("deprecated")`
#' `get_transform_information_from_stanfit()` was deprecated in MVBeliefUpdatr 0.1.0
#' and will be removed in 0.2.0. Please use [get_transform_information()] instead.
#' @seealso [get_transform_information()]
#' @keywords internal
#' @export
get_transform_information_from_stanfit <- function(...) {
  lifecycle::deprecate_warn(
    when = "0.1.0",
    what = "get_transform_information_from_stanfit()",
    with = "get_transform_information()"
  )
  get_transform_information(...)
}

#' Deprecated: get_transform_function_from_stanfit
#'
#' @description `r lifecycle::badge("deprecated")`
#' `get_transform_function_from_stanfit()` was deprecated in MVBeliefUpdatr 0.1.0
#' and will be removed in 0.2.0. Please use [get_transform_function()] instead.
#' @seealso [get_transform_function()]
#' @keywords internal
#' @export
get_transform_function_from_stanfit <- function(...) {
  lifecycle::deprecate_warn(
    when = "0.1.0",
    what = "get_transform_function_from_stanfit()",
    with = "get_transform_function()"
  )
  get_transform_function(...)
}

#' Deprecated: get_untransform_function_from_stanfit
#'
#' @description `r lifecycle::badge("deprecated")`
#' `get_untransform_function_from_stanfit()` was deprecated in MVBeliefUpdatr 0.1.0
#' and will be removed in 0.2.0. Please use [get_untransform_function()] instead.
#' @seealso [get_untransform_function()]
#' @keywords internal
#' @export
get_untransform_function_from_stanfit <- function(...) {
  lifecycle::deprecate_warn(
    when = "0.1.0",
    what = "get_untransform_function_from_stanfit()",
    with = "get_untransform_function()"
  )
  get_untransform_function(...)
}

#' Deprecated: get_staninput_from_stanfit
#'
#' @description `r lifecycle::badge("deprecated")`
#' `get_staninput_from_stanfit()` was deprecated in MVBeliefUpdatr 0.1.0
#' and will be removed in 0.2.0. Please use [get_staninput()] instead.
#' @seealso [get_staninput()]
#' @keywords internal
#' @export
get_staninput_from_stanfit <- function(...) {
  lifecycle::deprecate_warn(
    when = "0.1.0",
    what = "get_staninput_from_stanfit()",
    with = "get_staninput()"
  )
  get_staninput(...)
}

# get_exposure_category_statistic_from_stanfit <- get_exposure_category_statistic.ideal_adaptor_stanfit
# get_exposure_mean_from_stanfit <- get_exposure_category_mean.ideal_adaptor_stanfit
# get_exposure_css_from_stanfit <- get_exposure_category_css.ideal_adaptor_stanfit
# get_exposure_uss_from_stanfit <- get_exposure_category_uss.ideal_adaptor_stanfit
# get_exposure_cov_from_stanfit <- get_exposure_category_cov.ideal_adaptor_stanfit

#' Deprecated: get_test_data_from_stanfit
#'
#' @description `r lifecycle::badge("deprecated")`
#' `get_test_data_from_stanfit()` was deprecated in MVBeliefUpdatr 0.1.0
#' and will be removed in 0.2.0. Please use [get_test_data()] instead.
#' @seealso [get_test_data()]
#' @keywords internal
#' @export
get_test_data_from_stanfit <- function(...) {
  lifecycle::deprecate_warn(
    when = "0.1.0",
    what = "get_test_data_from_stanfit()",
    with = "get_test_data()"
  )
  get_test_data(...)
}

# get_original_variable_levels_from_stanfit <- get_staninput_variable_levels
# get_category_levels_from_stanfit <- get_category_levels
# get_group_levels_from_stanfit <- get_group_levels
# get_cue_levels_from_stanfit <- get_cue_levels
#
# get_expected_category_statistic_from_stanfit <- get_expected_category_statistic
# get_expected_mu_from_stanfit <- get_expected_mu
# get_expected_sigma_from_stanfit <- get_expected_sigma

#' Deprecated: add_ibbu_stanfit_draw
#'
#' @description `r lifecycle::badge("deprecated")`
#' `add_ibbu_stanfit_draw()` was deprecated in MVBeliefUpdatr 0.1.0
#' and will be removed in 0.2.0. Please use [get_draws()] instead.
#' @seealso [get_draws()]
#' @keywords internal
#' @export
add_ibbu_stanfit_draw <- function(...) {
  lifecycle::deprecate_warn(
    when = "0.1.0",
    what = "add_ibbu_stanfit_draw()",
    with = "get_draws()"
  )
  get_draws(...)
}

