#' Infer prior beliefs based on adaptation behavior
#'
#' This function takes exposure and test data, the names of the cue, category, and response columns (and
#' optionally group and/or block columns), and uses
#' stan to draw samples from the prior beliefs that match the behavior.
#'
#' @inheritParams compose_data_to_infer_prior_via_conjugate_ibbu_w_sufficient_stats
#' @param sample Should the model be fit and sampled from?
#' @param ... Additional parameters are passed to \code{\link[rstan]{sampling}}
#'
#' @return A \code{mvg_ibbu_stanfit} object with the fitted stan model.
#'
#' @seealso \code{\link{is.mvg_ibbu_stanfit}} for information about mvg_ibbu_stanfit objects,
#' \code{\link{add_ibbu_stanfit_draws}} to draw samples from the stanfit.
#' @export
infer_prior_beliefs <- function(
  exposure, test,
  cues, category, response, group, group.unique,
  center.observations = TRUE, scale.observations = TRUE, pca.observations = FALSE, pca.cutoff = 1,
  tau_scale = 10, L_omega_scale = 1,
  useMultivariateUpdating = if(length(cues) > 1) TRUE else FALSE,
  sample = TRUE, verbose = FALSE,
  ...) {
  data_list <- compose_data_to_infer_prior_via_conjugate_ibbu_w_sufficient_stats(
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
      tau_scale = tau_scale,
      L_omega_scale = L_omega_scale,
      verbose = verbose)

  message("Add compilation only options as in brms. and add in the data")

  if (sample) {
    if (useMultivariateUpdating) {
      fit <- rstan::stan(file = '../stancode/Updating_mvg_NIW_uninformativePriors_SufficientStats_wLapse',
                         data = data_list, ...)
    } else {
      fit <- rstan::sampling(beliefupdatr:::stanmodels[['conj_id_lapsing_sufficient_stats_fit']],
                             data = data_list, ...)
    }
  }

  fit %<>%
    as.mvg_ibbu_stanfit(data_list)

  return(fit)
}
