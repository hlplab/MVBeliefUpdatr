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
#'
#' Use \code{\link{infer_NIW_ideal_adaptor()}} instead, together with \code{\link{make_staninput_for_NIW_ideal_adaptor()}}.
#' @inheritParams make_staninput
#' @inheritParams fit_NIW_ideal_adaptor
#' @export
infer_prior_beliefs <- function(
  # arguments for make_staninput
  exposure, test,
  cues, category, response,
  group, group.unique = NULL,
  center.observations = TRUE, scale.observations = FALSE, pca.observations = FALSE, pca.cutoff = 1,
  lapse_rate = NULL, mu_0 = NULL, Sigma_0 = NULL,
  tau_scale = NULL, L_omega_eta = 1,
  split_loglik_per_observation = 0,
  stanmodel = NULL,
  verbose = FALSE,
  ...
) {
  # Currently the make_staninput function is creating both the transforms *and* the data.
  # That's a bit confusing and should probably be split up in the future into separate
  # functions.
  staninput <-
    make_staninput(
      exposure = exposure,
      test = test,
      cues = cues,
      category = category,
      response = response,
      group = group,
      group.unique = group.unique,
      center.observations = center.observations,
      scale.observations = scale.observations,
      pca.observations = pca.observations,
      pca.cutoff = pca.cutoff,
      lapse_rate = lapse_rate,
      mu_0 = mu_0,
      Sigma_0 = Sigma_0,
      tau_scale = tau_scale,
      L_omega_eta = L_omega_eta,
      split_loglik_per_observation = split_loglik_per_observation,
      use_univariate_updating = if (is.null(stanmodel)) { FALSE } else { stanmodel == 'uvg_conj_uninformative_priors_sufficient_stats_lapse'},
      verbose = verbose,
      model_type = "NIW_ideal_adaptor")

  fit_NIW_ideal_adaptor(staninput = staninput, stanmodel = stanmodel, ...)
}
