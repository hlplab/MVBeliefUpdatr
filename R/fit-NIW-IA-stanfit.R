#' Infer NIW ideal adaptor
#'
#' Infers a posterior distribution of NIW_ideal _adaptors from the input data using rstan/stan. The function can take
#' two types of inputs: an staninput list, as prepared by \code{\link{make_staninput}},
#' or the exposure and test data, the names of the cues, category, and response columns (and optionally group and/or block columns).
#'
#' @inheritParams make_staninput
#' @param staninput A list of the type that would be returned by \code{\link{make_staninput}}.
#'   This list can be provided *instead* of the arguments required by \code{make_staninput}.
#' @param untransform_fit Logical flag indicating whether the samples of the model should be transformed back
#'   into the original cue space by applying the untransform function. (default: \code{FALSE})
#' @param sample Should the model be fit and sampled from?
#' @param backend Character string naming the package to use as the backend for
#'   fitting the Stan model. Options are \code{"rstan"} (the default) or
#'   \code{"cmdstanr"}. Details on the
#'   \pkg{rstan} and \pkg{cmdstanr} packages are available at
#'   \url{https://mc-stan.org/rstan/} and \url{https://mc-stan.org/cmdstanr/},
#'   respectively.
#' @param file Either NULL or a character string. In the latter case, the fitted model object is saved
#'   via saveRDS in a file named after the string supplied in file. The .rds extension is added automatically.
#'   If the file already exists, the model from that file will be loaded and returned instead of refitting the model.
#'   As existing files won't be overwritten, you have to manually remove the file in order to refit and save the
#'   model under an existing file name. The file name is stored in the \code{NIW_ideal_adaptor_stanfit} object
#'   for later use. (default: \code{NULL})
#' @param file_compress Logical or a character string, specifying one of the
#'   compression algorithms supported by \code{\link{saveRDS}}. If the
#'   \code{file} argument is provided, this compression will be used when saving
#'   the fitted model object.
#' @param file_refit Modifies when the fit stored via the \code{file} argument
#'   is re-used. For \code{"never"} (default) the fit is always loaded if it
#'   exists and fitting is skipped. For \code{"always"} the model is always
#'   refitted. If set to \code{"on_change"}, model will be refit if data passed to
#'   Stan differ from
#'   what is stored in the file. Refit will not be triggered for changes in
#'   additional parameters of the fit (e.g., initial values, number of iterations,
#'   control arguments, ...).
#' @param silent Verbosity level between \code{0} and \code{2}.
#'   If \code{1} (the default), most of the
#'   informational messages of compiler and sampler are suppressed.
#'   If \code{2}, even more messages are suppressed. The actual
#'   sampling progress is still printed. Set \code{refresh = 0} to turn this off
#'   as well. If using \code{backend = "rstan"} you can also set
#'   \code{open_progress = FALSE} to prevent opening additional progress bars.
#' @param stan_model_args A \code{list} of further arguments passed to
#'   \code{\link[rstan:stan_model]{rstan::stan_model}} for \code{backend =
#'   "rstan"} or to \code{cmdstanr::cmdstan_model} for \code{backend =
#'   "cmdstanr"}, which allows to change how models are compiled.
#' @param stanmodel Name of stanmodel that should be used. Overrides any default selection.
#' @param ... Additional parameters are passed to \code{\link[rstan]{sampling}}
#'
#' @return An object of class \code{NIW_ideal_adaptor_stanfit} with the fitted stan model. In interpreting the
#'   inferred parameters, it should be kept in mind that the \emph{inferred} scatter matrix S_0 includes
#'   variability from internal perceptual and/or external environmental noise, \emph{in addition} to the motor
#'   noise that is reflected in production data. This also implies that, \strong{if \code{Sigma_0} is provided
#'   by the user it should be convolved with perceptual noise}. This is particularly important if the data you're
#'   fitting contains test phases without exposure (e.g., pre-exposure tests). Make sure to read the notes about
#'   the \code{Sigma_0} argument in the help page on \code{\link{make_staninput}}. Use
#'   \code{methods(class = "NIW_ideal_adaptor_stanfit")} for an overview on available methods.
#'
#' @seealso \code{\link{is.NIW_ideal_adaptor_stanfit}} for information about NIW_ideal_adaptor_stanfit objects,
#' \code{\link{get_draws}} to draw samples from the stanfit.
#'
#' @importFrom rstan nlist
#' @rdname infer_NIW_ideal_adaptor
#' @export
infer_NIW_ideal_adaptor <- function(
  # arguments for make_staninput
  exposure = NULL, test = NULL,
  cues = NULL, category = NULL, response = NULL, group = NULL, group.unique = NULL,
  center.observations = TRUE, scale.observations = FALSE, pca.observations = FALSE, pca.cutoff = 1,
  lapse_rate = NULL, mu_0 = NULL, Sigma_0 = NULL,
  tau_scale = NULL, L_omega_eta = 1,
  split_loglik_per_observation = 0,
  # additional arguments:
  staninput = NULL,
  # set to default TRUE (and change documentation accordingly) once implemented:
  untransform_fit = FALSE,
  sample = TRUE,
  file = NULL, file_refit = "never", file_compress = T,
  # included for later use
  stanvars = NULL, backend = "rstan", save_pars = NULL, basis = NULL,
  init = NULL, control = NULL, silent = 1, stan_model_args = list(),
  # Stuff to be deprecated in the future
  stanmodel = NULL,
  verbose = FALSE,
  ...
) {
  # optionally load NIW_ideal_adaptor_stanfit from file
  # Loading here only when we should directly load the file.
  # The "on_change" option needs more information
  file_refit <- match.arg(file_refit, file_refit_options())
  if (!is.null(file) && file_refit == "never") {
    fit <- read_NIW_ideal_adaptor_stanfit(file)
    if (!is.null(fit)) {
      return(fit)
    }
  }

  assert_that(
    all(any(all(!is.null(exposure), !is.null(test)), !is.null(staninput)), !all(!is.null(exposure), !is.null(test), !is.null(staninput))),
    msg = "You can either specify an staninput argument or exposure and test data, but not both.")

  # Some initial checking of the staninput list information. This should probably be expanded and then collected into a function
  # that checks what structure compose_ returns.
  if (!is.null(staninput)) {
    assert_that(all(!is.null(staninput$staninput), !is.null(staninput$transform_information)),
                msg = "The staninput list must contain both staninput and transform_information.")
    assert_that(is.list(staninput$transform_information))
    assert_that(all(c("transform.function", "untransform.function") %in% names(staninput$transform_information)),
                msg = "If not NULL, transform_information must be a list that contains both a transform and an untransform function.")
  }

  if (!is.null(stanmodel)) {
    assert_that(!is.null(stanmodels[[stanmodel]]),
                msg = paste("The specified stanmodel does not exist. Allowable models include:", names(MVBeliefUpdatr:::stanmodels)))
  }

  if (is.null(staninput)) {
    staninput <-
      # Currently the make_staninput function is creating both the transforms *and* the data.
      # That's a bit confusing and should probably be split up in the future into separate
      # functions.
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
  }

  data <-
    bind_rows(
      exposure[, c(group.unique, group, category, cues)] %>%
        drop_na() %>%
        mutate(Phase = "exposure"),
      exposure[, c(group.unique, group, response, cues)] %>%
        drop_na() %>%
        mutate(Phase = "test")) %>%
    relocate(Phase, all_of(c(group.unique, group, category, cues, response)))

  # Check whether model actually needs to be refit
  if (!is.null(file) && file_refit == "on_change") {
    x_from_file <- read_NIW_ideal_adaptor_stanfit(file)
    if (!is.null(x_from_file)) {
      needs_refit <-
        stanfit_needs_refit(
          x_from_file,
          current_version = get_current_versions(),
          data = data, staninput = staninput,
          silent = silent)
      if (!needs_refit) {
        return(x_from_file)
      }
    }
  }


  # Clean-up x_mean and x_ss for groups without exposure data. For reasons laid out in
  # get_sufficient_statistics_as_list_of_arrays, we had to set these means and sums of
  # squares to arbitrary values (since Stan doesn't accept typed NAs). But this can
  # create confusion when users try to retrieve the exposure statistics for those groups.
  # Here we're thus setting them to NAs.
  staninput$staninput$x_mean[staninput$staninput$N == 0] <- NA
  staninput$staninput$x_ss[staninput$staninput$N == 0] <- NA
  fit <-
    NIW_ideal_adaptor_stanfit(
      data = data,
      staninput = staninput$staninput,
      stanvars = stanvars,
      save_pars = save_pars,
      backend = backend,
      stan_args = nlist(init, silent, control, stan_model_args),
      transform_information = staninput$transform_information,
      basis = basis,
      file = file)

  current_default_modelname <- 'mvg_conj_sufficient_stats_lapse_cholesky'
  stanfit <- NULL
  if (sample) {
    if (is.null(stanmodel) || stanmodel == current_default_modelname) {
      # Parameters *not* to store
      exclude_pars <-
        c("lapse_rate_param",
          "m_0_param", "m_0_tau", "m_0_tau_param", "m_0_L_omega", "m_0_L_omega_param",
          "tau_0_param", "L_omega_0_param", "L_S_0", "L_S_n", "L_t_scale",
          "p_test_conj", "log_p_test_conj")

      stanfit <-
        sampling(
          MVBeliefUpdatr:::stanmodels[[current_default_modelname]],
          data = staninput$staninput,
          pars = exclude_pars, include = FALSE,
          show_messages = !silent,
          ...)
    } else if (stanmodel %in% names(MVBeliefUpdatr:::stanmodels)) {
        stanfit <-
          sampling(
            MVBeliefUpdatr:::stanmodels[[stanmodel]],
            data = staninput$staninput,
            show_messages = !silent,
            ...)

    }

    if (is.null(stanfit)) {
      stop("Sampling failed.")
    } else {
      fit %<>% set_stanfit(stanfit)
      fit %<>% recover_types.NIW_ideal_adaptor_stanfit()
    }
  }

  if (untransform_fit) {
    message("untransform_fit has not yet been implemented. Returning unchanged fit.")
  }

  if (!is.null(fit) & !is.null(file)) {
    fit <- write_NIW_ideal_adaptor_stanfit(fit, file, compress = file_compress)
  }

  return(fit)
}

