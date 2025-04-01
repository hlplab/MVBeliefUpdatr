#' DEPRECATED: get_transform_information_from_stanfit
#' @export
get_transform_information_from_stanfit <- function(...) get_transform_information.NIW_ideal_adaptor_stanfit(...)

#' DEPRECATED: get_transform_function_from_stanfit
#' @export
get_transform_function_from_stanfit <- function(...) get_transform_function.NIW_ideal_adaptor_stanfit(...)

#' DEPRECATED: get_untransform_function_from_stanfit
#' @export
get_untransform_function_from_stanfit <- function(...) get_untransform_function.NIW_ideal_adaptor_stanfit(...)

#' DEPRECATED: get_staninput_from_stanfit
#' @export
get_staninput_from_stanfit <- function(...) get_staninput.NIW_ideal_adaptor_stanfit(...)

# get_exposure_category_statistic_from_stanfit <- get_exposure_category_statistic.NIW_ideal_adaptor_stanfit
# get_exposure_mean_from_stanfit <- get_exposure_category_mean.NIW_ideal_adaptor_stanfit
# get_exposure_css_from_stanfit <- get_exposure_category_css.NIW_ideal_adaptor_stanfit
# get_exposure_uss_from_stanfit <- get_exposure_category_uss.NIW_ideal_adaptor_stanfit
# get_exposure_cov_from_stanfit <- get_exposure_category_cov.NIW_ideal_adaptor_stanfit

#' DEPRECATED: get_test_data_from_stanfit
#' @export
get_test_data_from_stanfit <- function(...) get_test_data(...)

# get_original_variable_levels_from_stanfit <- get_staninput_variable_levels
# get_category_levels_from_stanfit <- get_category_levels
# get_group_levels_from_stanfit <- get_group_levels
# get_cue_levels_from_stanfit <- get_cue_levels
#
# get_expected_category_statistic_from_stanfit <- get_expected_category_statistic
# get_expected_mu_from_stanfit <- get_expected_mu
# get_expected_sigma_from_stanfit <- get_expected_sigma

#' DEPRECATED: add_ibbu_stanfit_draw
#' @export
add_ibbu_stanfit_draw <- function(...) get_draws(...)

#' DEPRECATED: Infer prior beliefs
#' @export
infer_prior_beliefs <- function(...) infer_NIW_ideal_adaptor(...)
