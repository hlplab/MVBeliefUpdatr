#' Infer prior beliefs based on adaptation behavior
#'
#' This function takes exposure and test data, the names of the cues, category, and response columns (and
#' optionally group and/or block columns), and uses
#' stan to draw samples from the prior beliefs that match the behavior.
#'
#' @inheritParams compose_data_to_infer_prior_via_conjugate_ibbu_w_sufficient_stats
#' @param sample Should the model be fit and sampled from?
#' @param file Either NULL or a character string. In the latter case, the fitted model object is saved
#' via saveRDS in a file named after the string supplied in file. The .rds extension is added automatically.
#' If the file already exists, the model from that file will be loaded and returned instead of refitting the model.
#' As existing files won't be overwritten, you have to manually remove the file in order to refit and save the
#' model under an existing file name. The file name is stored in the NIW_ibbu_stanfit object for later usage.
#' (default `NULL`)
#' @param use_multivariate_updating Should multivariate updating be used? By default this option will be
#' selecting if the relevant cues (after transformations, including PCA, if selected) have more than 1
#' dimension.
#' @param model Name of stanmodel that should be used. Overrides any default selection.
#' @param ... Additional parameters are passed to \code{\link[rstan]{sampling}}
#'
#' @return \code{NIW_ibbu_stanfit} object with the fitted stan model.
#'
#' @seealso \code{\link{is.NIW_ibbu_stanfit}} for information about NIW_ibbu_stanfit objects,
#' \code{\link{add_ibbu_stanfit_draws}} to draw samples from the stanfit.
#' @export
infer_prior_beliefs <- function(
  exposure, test, data_list = NULL,
  cues, category, response, group, group.unique,
  center.observations = TRUE, scale.observations = TRUE, pca.observations = FALSE, pca.cutoff = 1,
  m_0 = NULL, S_0 = NULL,
  tau_scale = 0, L_omega_scale = 0,
  sample = TRUE, file = NULL, model = NULL, verbose = FALSE,
  use_multivariate_updating = NULL,
  ...) {
  if (!is.null(file)) {
    x <- read_NIW_ibbu_stanfit(file)
    if (!is.null(x)) {
      return(x)
    }
  }
  if (!is.null(model)) assert_that(!is.null(stanmodels[[model]]),
                                   msg = paste("The specified stanmodel does not exist. Allowable models include:",
                                               names(MVBeliefUpdatr:::stanmodels)))
  if (!is.null(data_list)) {
      transform <- NULL
  } else {
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
      m_0 = m_0,
      S_0 = S_0,
      tau_scale = tau_scale,
      L_omega_scale = L_omega_scale,
      verbose = verbose)

    transform <-
      transform_cues(
        data = exposure,
        cues = cues,
        center = center.observations,
        scale = scale.observations,
        pca = pca.observations,
        return.transformed.data = F,
        return.transform.parameters = F,
        return.transform.function = T,
      return.untransform.function = T)
  }

  if (is.null(use_multivariate_updating))
    use_multivariate_updating = if (is.null(data_list$K)) FALSE else TRUE
  message("change to keep samples only from relevant parameters.")

  if (sample) {
    if (!is.null(model)) {
      fit <- sampling(MVBeliefUpdatr:::stanmodels[[model]],
                      data = data_list, ...)
    } else if (use_multivariate_updating) {
      fit <- sampling(MVBeliefUpdatr:::stanmodels[['mvg_conj_uninformative_priors_sufficient_stats_lapse']],
                         data = data_list, ...)
    } else {
      message("There might be an issue with the compose_data function for univariate models. look into it.")
      fit <- sampling(MVBeliefUpdatr:::stanmodels[['uvg_conj_uninformative_priors_sufficient_stats_lapse']],
                             data = data_list, ...)
    }

  if (is.null(fit)) stop("Sampling failed.")
  fit %<>% as.NIW_ibbu_stanfit(data_list, transform)
  fit %<>%
    recover_types(
      crossing(
        category = factor(colnames(data_list$z_test_counts), levels = colnames(data_list$z_test_counts)),
        group = factor(attr(data_list$y_test, "levels"), levels = attr(data_list$y_test, "levels")),
        cue = factor(names(data_list$x_test), levels = names(data_list$x_test)),
        cue2 = cue
      )
    )
  }

  if (!is.null(file)) {
    write_NIW_ibbu_stanfit(fit, file)
  }

  return(fit)
}
