#' Infer NIW ideal adaptor
#'
#' Infers a posterior distribution of NIW_ideal _adaptors from the input data using rstan/stan. The function can take
#' two types of inputs: an input list, as prepared by \code{\link{compose_data_to_infer_NIW_ideal_adaptor}},
#' or the exposure and test data, the names of the cues, category, and response columns (and optionally group and/or block columns).
#'
#' @inheritParams compose_data_to_infer_NIW_ideal_adaptor
#' @param untransform_fit Logical flag indicating whether the samples of the model should be transformed back
#' into the original cue space by applying the untransform function. (default: \code{FALSE})
#' @param input A list of the type that would be returned by \code{\link{compose_data_to_infer_NIW_ideal_adaptor}}.
#' This list can be provided *instead* of the arguments required by \code{compose_data_to_infer_NIW_ideal_adaptor}.
#' @param sample Should the model be fit and sampled from?
#' @param file Either NULL or a character string. In the latter case, the fitted model object is saved
#' via saveRDS in a file named after the string supplied in file. The .rds extension is added automatically.
#' If the file already exists, the model from that file will be loaded and returned instead of refitting the model.
#' As existing files won't be overwritten, you have to manually remove the file in order to refit and save the
#' model under an existing file name. The file name is stored in the \code{NIW_ideal_adaptor_stanfit} object
#' for later use. (default: \code{NULL})
#' @param model Name of stanmodel that should be used. Overrides any default selection.
#' @param use_univariate_updating Should legacy univariate updating be used? Will throw an error if used in
#' conjunction with multiple cues. (default: \code{FALSE})
#' @param ... Additional parameters are passed to \code{\link[rstan]{sampling}}
#'
#' @return \code{NIW_ideal_adaptor_stanfit} object with the fitted stan model. In interpreting the inferred parameters, it should
#' be kept in mind that the \emph{inferred} scatter matrix S_0 includes variability from internal perceptual and/or
#' external environmental noise, \emph{in addition} to the motor noise that is reflected in production data. This also
#' implies that, \strong{if \code{Sigma_0} is provided by the user it should be convolved with perceptual noise}. This is particularly important
#' if the data you're fitting contains test phases without exposure (e.g., pre-exposure tests). Make sure to read the notes about the
#' \code{Sigma_0} argument in the help page on \code{\link{compose_data_to_infer_NIW_ideal_adaptor}}.
#'
#' @seealso \code{\link{is.NIW_ideal_adaptor_stanfit}} for information about NIW_ideal_adaptor_stanfit objects,
#' \code{\link{get_draws}} to draw samples from the stanfit.
#' @rdname infer_NIW_ideal_adaptor
#' @export
infer_NIW_ideal_adaptor <- function(
  exposure = NULL, test = NULL,
  cues = NULL,  category = NULL, response = NULL,
  group = NULL, group.unique = NULL,
  center.observations = TRUE, scale.observations = TRUE, pca.observations = FALSE, pca.cutoff = 1,
  lapse_rate = NULL, mu_0 = NULL, Sigma_0 = NULL,
  tau_scale = 0, L_omega_scale = 0,
  split_loglik_per_observation = 0,
  # set to default TRUE (and change documentation accordingly) once implemented:
  untransform_fit = FALSE,
  input = NULL,
  sample = TRUE,
  file = NULL,
  model = NULL,
  use_univariate_updating = FALSE,
  verbose = FALSE,
  ...) {
  assert_that(
    all(any(all(!is.null(exposure), !is.null(test)), !is.null(input)), !all(!is.null(exposure), !is.null(test), !is.null(input))),
    msg = "You can either specify an input argument or exposure and test data, but not both.")

  # Some initial checking of the input list information. This should probably be expanded and then collected into a function
  # that checks what structure compose_ returns.
  if (!is.null(input)) {
    assert_that(all(!is.null(input$data_list), !is.null(input$transform_information)),
                msg = "The input list must contain both data_list and transform_information.")
    assert_that(is.list(input$transform_information))
    assert_that(all(c("transform.function", "untransform.function") %in% names(input$transform_information)),
                msg = "If not NULL, transform_information must be a list that contains both a transform and an untransform function.")
  }

  assert_that(is.logical(use_univariate_updating))
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

  if (is.null(input)) {
    input <-
      compose_data_to_infer_NIW_ideal_adaptor(
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
        split_loglik_per_observation = split_loglik_per_observation,
        use_univariate_updating = use_univariate_updating,
        verbose = verbose)
  }

  if (use_univariate_updating)
    assert_that(
      is.null(input$data_list$K),
      msg = "Univariate updating cannot be used with multiple cues.")

  if (sample) {
    # Parameters *not* to store
    pars <- c("m_0_tau", "m_0_L_omega", "tau_0_param", "t_scale", "p_test_conj", "log_p_test_conj")

    if (!is.null(model)) {
      fit <- sampling(MVBeliefUpdatr:::stanmodels[[model]],
                      data = input$data_list, ...)
    } else if (use_univariate_updating) {
      message("There might be an issue with the compose_data function for univariate models. look into it.")
      fit <- sampling(MVBeliefUpdatr:::stanmodels[['uvg_conj_uninformative_priors_sufficient_stats_lapse']],
                      data = input$data_list, ...)
    } else {
      fit <- sampling(MVBeliefUpdatr:::stanmodels[['mvg_conj_sufficient_stats_lapse']],
                      data = input$data_list, pars = pars, include = F, ...)
    }

    if (is.null(fit)) stop("Sampling failed.")
    # Clean-up x_mean and x_ss for groups without exposure data. For reasons laid out in
    # get_sufficient_statistics_as_list_of_arrays, we had to set these means and sums of
    # squares to arbitrary values (since Stan doesn't accept typed NAs). But this can
    # create confusion when users try to retrieve the exposure statistics for those groups.
    # Here we're thus setting them to NAs.
    input$data_list$x_mean[input$data_list$N == 0] <- NA
    input$data_list$x_ss[input$data_list$N == 0] <- NA
    fit %<>% as.NIW_ideal_adaptor_stanfit(input_data = input$data_list, transform_information = input$transform_information)
  } else fit <- NULL

  if (untransform_fit) {
    message("untransform_fit has not yet been implemented. Returning unchanged fit.")
  }

  if (!is.null(fit) & !is.null(file)) {
    write_NIW_ideal_adaptor_stanfit(fit, file)
  }

  return(fit)
}

infer_prior_beliefs <- infer_NIW_ideal_adaptor
