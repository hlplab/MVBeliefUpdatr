#' Infer prior beliefs based on adaptation behavior
#'
#' This function takes exposure and test data, the names of the cues, category, and response columns (and
#' optionally group and/or block columns), and uses
#' stan to draw samples from the prior beliefs that match the behavior.
#'
#' @inheritParams compose_data_to_infer_prior_via_conjugate_ibbu_w_sufficient_stats
#' @param data_list A list of the type that would be returned by \code{\link[compose_data]{compose_data_to_infer_prior_via_conjugate_ibbu_w_sufficient_stats}}.
#' This list can be provided instead of the individual inputs that are handed to \code{compose_data_to_infer_prior_via_conjugate_ibbu_w_sufficient_stats}.
#' @param transform_information Optionally, a list of transform paramters, transform function, and corresponding untransform
#' function of the type returned by \code{\link[transform_cues]{transform_cues}}. These functions
#' will be attached to the \code{NIW_ideal_adaptor_stanfit} object. This object will override any transform
#' information that is created by \code{\link[compose_data]{compose_data_to_infer_prior_via_conjugate_ibbu_w_sufficient_stats}}.
#' @param untransform_fit Logical flag indicating whether the MCMC samples of the model should be transformed back
#' into the original cue space by applying the untransform function. (default: `TRUE`)
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
  m_0 = NULL, S_0 = NULL,
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
        m_0 = m_0,
        S_0 = S_0,
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
  message("message to developer: change to keep samples only from relevant parameters.")

  if (sample) {
    pars <- c()
    if (!is.null(m_0)) {
      pars <- c(pars)

      # m_0_tau[1]              8.5e+307     NaN     Inf  9.0e+306  4.4e+307  7.8e+307  1.2e+308  1.8e+308   NaN  NaN
      # m_0_tau[2]              1.0e+308     NaN     Inf  9.8e+306  6.5e+307  1.1e+308  1.5e+308  1.8e+308   NaN  NaN
      # m_0_L_omega[1,1]         1.0e+00     NaN 0.0e+00   1.0e+00   1.0e+00   1.0e+00   1.0e+00   1.0e+00   NaN  NaN
      # m_0_L_omega[1,2]         0.0e+00     NaN 0.0e+00   0.0e+00   0.0e+00   0.0e+00   0.0e+00   0.0e+00   NaN  NaN
      # m_0_L_omega[2,1]        -5.0e-02 9.0e-02 2.0e-01  -3.0e-01  -1.9e-01  -8.0e-02   3.0e-02   5.9e-01     4  1.9
      # m_0_L_omega[2,2]         9.8e-01 1.0e-02 4.0e-02   8.1e-01   9.8e-01   9.9e-01   1.0e+00   1.0e+00    13  1.3
      # tau_0_param[1,1]         7.2e+01 4.8e+01 6.8e+01   3.0e+01   3.2e+01   3.4e+01   7.2e+01   2.0e+02     2 24.5
      # tau_0_param[1,2]         4.7e+01 1.8e+01 2.6e+01   2.8e+01   3.2e+01   3.4e+01   4.9e+01   9.6e+01     2 14.6
      # tau_0_param[2,1]         2.7e+01 5.7e+00 8.1e+00   1.3e+01   2.4e+01   3.0e+01   3.3e+01   3.6e+01     2  7.2
      # tau_0_param[2,2]         6.6e+01 2.0e+01 2.9e+01   1.6e+01   4.9e+01   7.7e+01   8.9e+01   9.6e+01     2 12.9
      # m_0[1,1]                 8.9e-01 0.0e+00 0.0e+00   8.9e-01   8.9e-01   8.9e-01   8.9e-01   8.9e-01     2  1.0
      # m_0[1,2]                 9.1e+00 0.0e+00 0.0e+00   9.1e+00   9.1e+00   9.1e+00   9.1e+00   9.1e+00     2  1.0
      # m_0[2,1]                 3.4e+00 0.0e+00 0.0e+00   3.4e+00   3.4e+00   3.4e+00   3.4e+00   3.4e+00     2  1.0
      # m_0[2,2]                 9.3e+00     NaN 0.0e+00   9.3e+00   9.3e+00   9.3e+00   9.3e+00   9.3e+00   NaN  1.0


    }

    if (!is.null(model)) {
      fit <- sampling(MVBeliefUpdatr:::stanmodels[[model]],
                      data = data_list, ...)
    } else if (use_multivariate_updating) {
      fit <- sampling(MVBeliefUpdatr:::stanmodels[['mvg_conj_sufficient_stats_lapse']],
                         data = data_list, ...)
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
