#' Infer prior beliefs based on adaptation behavior
#'
#' This function takes exposure and test data, the names of the cues, category, and response columns (and
#' optionally group and/or block columns), and uses
#' stan to draw samples from the prior beliefs that match the behavior.
#'
#' @inheritParams compose_data_to_infer_prior_via_conjugate_ibbu_w_sufficient_stats
#' @param data_list A list of the type that would be returned by \code{\link[compose_data]{compose_data_to_infer_prior_via_conjugate_ibbu_w_sufficient_stats}}.
#' This list can be provided instead of the individual inputs that are handed to \code{compose_data_to_infer_prior_via_conjugate_ibbu_w_sufficient_stats}.
#' @param untransform_fit Logical flag indicating whether the samples of the model should be transformed back
#' into the original cue space by applying the untransform function. (default: `TRUE`)
#' @param transform_information Optionally, a list of transform paramters, transform function, and corresponding untransform
#' function of the type returned by \code{\link[transform_cues]{transform_cues}}. These functions
#' will be attached to the \code{NIW_ideal_adaptor_stanfit} object. This object will override any transform
#' information that is created by \code{\link[compose_data]{compose_data_to_infer_prior_via_conjugate_ibbu_w_sufficient_stats}}.
#' @param sample Should the model be fit and sampled from?
#' @param file Either NULL or a character string. In the latter case, the fitted model object is saved
#' via saveRDS in a file named after the string supplied in file. The .rds extension is added automatically.
#' If the file already exists, the model from that file will be loaded and returned instead of refitting the model.
#' As existing files won't be overwritten, you have to manually remove the file in order to refit and save the
#' model under an existing file name. The file name is stored in the \code{NIW_ideal_adaptor_stanfit} object
#' for later use. (default `NULL`)
#' @param model Name of stanmodel that should be used. Overrides any default selection.
#' @param use_multivariate_updating Should multivariate updating be used? By default this option will be
#' selecting if the relevant cues (after transformations, including PCA, if selected) have more than 1
#' dimension.
#' @param ... Additional parameters are passed to \code{\link[rstan]{sampling}}
#'
#' @return \code{NIW_ideal_adaptor_stanfit} object with the fitted stan model.
#'
#' @seealso \code{\link{is.NIW_ideal_adaptor_stanfit}} for information about NIW_ideal_adaptor_stanfit objects,
#' \code{\link{add_ibbu_stanfit_draws}} to draw samples from the stanfit.
#' @export
infer_prior_beliefs <- function(
  exposure, test,
  cues, category, response, group, group.unique = NULL,
  center.observations = TRUE, scale.observations = TRUE, pca.observations = FALSE, pca.cutoff = 1, untransform_fit = TRUE,
  lapse_rate = NULL, mu_0 = NULL, Sigma_0 = NULL,
  tau_scale = 0, L_omega_scale = 0,
  data_list = NULL, transform_information = NULL,
  sample = TRUE, file = NULL, model = NULL, use_multivariate_updating = NULL,
  verbose = FALSE,
  ...) {
  if (!is.null(file)) {
    fit <- read_NIW_ideal_adaptor_stanfit(file)
    if (!is.null(fit)) {
      return(fit)
    }
  }
  if (!is.null(model)) {
    assert_that(!is.null(stanmodels[[model]]),
                msg = paste("The specified stanmodel does not exist. Allowable models include:", names(MVBeliefUpdatr:::stanmodels)))
  }
  if (!is.null(transform_information)) {
    assert_that(is.list(transform_information))
    assert_that(all(c("transform.function", "untransform.function") %in% names(transform_information)),
                msg = "If not NULL, transform_information must be a list that contains both a transform and an untransform function.")
  }

  if (is.null(data_list)) {
    data_list <-
      compose_data_to_infer_prior_via_conjugate_ibbu_w_sufficient_stats(
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
        L_omega_scale = L_omega_scale,
        verbose = verbose)

    if (is.null(transform_information)) {
      transform_information <-
        transform_cues(
          data = exposure,
          cues = cues,
          center = center.observations,
          scale = scale.observations,
          pca = pca.observations,
          return.transformed.data = F,
          return.transform.parameters = T,
          return.transform.function = T,
          return.untransform.function = T)
    }
  }

  if (is.null(use_multivariate_updating))
    use_multivariate_updating = if (is.null(data_list$K)) FALSE else TRUE

  if (sample) {
    # Parameters *not* to store
    pars <- c("m_0_tau", "m_0_L_omega", "tau_0_param", "t_scale", "p_test_conj", "log_p_test_conj")

    if (!is.null(model)) {
      fit <- sampling(MVBeliefUpdatr:::stanmodels[[model]],
                      data = data_list, ...)
    } else if (use_multivariate_updating) {
      fit <- sampling(MVBeliefUpdatr:::stanmodels[['mvg_conj_sufficient_stats_lapse']],
                         data = data_list, pars = pars, include = F, ...)
    } else {
      message("There might be an issue with the compose_data function for univariate models. look into it.")
      fit <- sampling(MVBeliefUpdatr:::stanmodels[['uvg_conj_uninformative_priors_sufficient_stats_lapse']],
                             data = data_list, ...)
    }

    if (is.null(fit)) stop("Sampling failed.")
    fit %<>% as.NIW_ideal_adaptor_stanfit(input_data = data_list, transform_information = transform_information)
  } else fit <- NULL

  if (untransform_fit) {
    message("untransform_fit has not yet been implemented. Returning unchanged fit.")
  }

  if (!is.null(fit) & !is.null(file)) {
    write_NIW_ideal_adaptor_stanfit(fit, file)
  }

  return(fit)
}
