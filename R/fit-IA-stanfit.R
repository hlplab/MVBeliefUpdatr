#' Fit ideal adaptor
#'
#' Infers a prior and posterior distribution of ideal adaptors from the input data using Stan. Currently, three
#' types of ideal adaptor models are available, each using a conjugate prior over Gaussian or multivariate
#' Gaussian categories:
#'
#' \itemize{
#'  \item{\code{NIW_ideal_adaptor}:} A Normal-Inverse-Wishart (NIW) prior over the ideal adaptor. This is the
#'    default model with multivariate Gaussian categories. Accepts univariate and multivariate input, though
#'    the NIX model should be faster for univariate input.
#'  \item{\code{NIX_ideal_adaptor}:} A Normal-Inverse-Chisquare (NIX) prior over the ideal adaptor with with
#'    univariate Gaussian categories. Accepts only univariate input.
#'  \item{\code{MNIX_ideal_adaptor}:} Separate NIXs for each of multiple cues that are integrated over during
#'    categorization(cue integration) assuming ideal cue weights based on the relative informativity of each cue.
#'    Accepts univariate and multivariate input, though the NIX model should be faster for univariate input.
#' }
#'
#' @param staninput A list of the type returned by \code{\link{make_staninput}}.
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
#'   model under an existing file name. The file name is stored in the \code{ideal_adaptor_stanfit} object
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
#'   \code{\link[rstan]{stan_model}} for \code{backend =
#'   "rstan"} or to \code{cmdstanr::cmdstan_model} for \code{backend =
#'   "cmdstanr"}, which allows to change how models are compiled.
#' @param stanmodel Name of stanmodel that should be used. Overrides any default selection.
#' @param rename For internal use only.
#' @param ... Additional parameters are passed to \code{\link[rstan]{sampling}}.
#'
#' @return An object of class \code{ideal_adaptor_stanfit} with the fitted stan model.
#'
#' @details
#'   In interpreting the inferred parameters, it should be kept in mind that the \emph{inferred} scatter matrix
#'   (e.g., S_0) for the NIW model) includes
#'   variability from internal perceptual and/or external environmental noise, \emph{in addition} to the motor
#'   noise that is reflected in production data. This also implies that, \strong{if \code{Sigma_0} is provided
#'   by the user it should be arguably convolved with an estimate of perceptual noise}. This is particularly
#'   important if the data you're fitting contains test phases without exposure (e.g., pre-exposure tests).
#'
#'   Make sure to read the notes about the \code{Sigma_0} argument in the help page on \code{\link{make_staninput}}.
#'   Use \code{methods(class = "ideal_adaptor_stanfit")} for an overview on available methods.
#'
#' @seealso \code{\link{is.ideal_adaptor_stanfit}} for information about ideal_adaptor_stanfit objects,
#' \code{\link{get_draws}} to draw samples from the stanfit.
#'
#' @importFrom rstan nlist
#' @export
fit_ideal_adaptor <- function(
  staninput,
  file = NULL, file_refit = "never", file_compress = T,
  # included for later use
  stanvars = NULL, backend = "rstan", save_pars = NULL, basis = NULL,
  chains = 4, iter = 2000, warmup = 1000,
  init = "random", control = NULL,
  silent = 1, verbose = F,
  stan_model_args = list(),
  stanmodel = NULL,
  # Stuff to be deprecated in the future
  rename = T,
  ...
) {
  # optionally load ideal_adaptor_stanfit from file
  # Loading here only when we should directly load the file.
  # The "on_change" option needs more information
  file_refit <- match.arg(file_refit, file_refit_options())
  if (!is.null(file) && file_refit == "never") {
    fit <- read_ideal_adaptor_stanfit(file)
    if (!is.null(fit)) {
      if (silent == 0) message("Loading existing model from file.")
      return(fit)
    }
  }

  # extract information from staninput
  assert_that(is.ideal_adaptor_staninput(staninput, verbose = verbose))
  data <- staninput$data
  transform_information <- staninput$transform_information
  staninput <- staninput$staninput

  if (!is.null(stanmodel)) {
    assert_that(!is.null(stanmodels[[stanmodel]]),
                msg = paste("The specified stanmodel does not exist. Allowable models include:", paste(names(MVBeliefUpdatr:::stanmodels), collapse = ", ")))
  }

  # Check whether model actually needs to be refit
  if (!is.null(file) && file_refit == "on_change") {
    x_from_file <- read_ideal_adaptor_stanfit(file)
    if (!is.null(x_from_file)) {
      needs_refit <-
        stanfit_needs_refit(
          x_from_file,
          current_version = get_current_versions(),
          data = data, staninput = staninput,
          silent = silent, verbose = verbose)
      if (!needs_refit) {
        if (silent == 0) message("No refitting needed. Loading existing model from file.")
        return(x_from_file)
      }
    }
  }

  fit <-
    ideal_adaptor_stanfit(
      data = data,
      staninput = staninput,
      stanvars = stanvars,
      save_pars = save_pars,
      backend = backend,
      stan_args = nlist(init, silent, control, stan_model_args),
      transform_information = transform_information,
      basis = basis,
      file = file)

  # Check that staninput has at least two categories (fitting with one category makes no sense)
  if (get_staninput(fit, which = "transformed")$M < 2) stop("staninput must have at least two categories.")

  stanfit <- NULL
  if (chains > 0 & iter > 0) {
    # Parameters *not* to store (used with include = F below)
    # Future file reduction could be achieved via the shredder package for stanfit
    # post-processing (https://github.com/yonicd/shredder)
    exclude_pars <-
      c("lapse_rate_param",
        "m_0_param", "m_0_tau", "m_0_tau_param", "m_0_L_omega", "m_0_L_omega_param",
        "tau_0_param", "L_omega_0_param", "L_S_0", "L_S_n", "L_t_scale",
        "p_test_conj", "log_p_test_conj")

    if (is.null(stanmodel)) {
      current_default_modelname <- 'NIW_ideal_adaptor'
      stanfit <-
        sampling(
          MVBeliefUpdatr:::stanmodels[[current_default_modelname]],
          data = get_staninput(fit, which = "transformed"),
          check_data = TRUE,
          pars = exclude_pars, include = FALSE,
          chains = chains, iter = iter, warmup = warmup,
          init = init, control = control,
          show_messages = !silent,
          ...)
    } else if (stanmodel %in% names(MVBeliefUpdatr:::stanmodels)) {
        stanfit <-
          sampling(
            MVBeliefUpdatr:::stanmodels[[stanmodel]],
            data = get_staninput(fit, which = "transformed"),
            check_data = TRUE,
            pars = exclude_pars, include = FALSE,
            chains = chains, iter = iter, warmup = warmup,
            init = init, control = control,
            show_messages = !silent,
            ...)

    }

    if (!contains_draws(stanfit)) {
      stop2("Sampling failed.")
    } else {
      fit %<>% set_stanfit(stanfit)
      fit %<>% recover_types.ideal_adaptor_stanfit()

      if (rename) fit %<>% rename_pars()
    }
  } else if (!silent) message("No sampling requested. Returning empty model object.")

  if (!is.null(fit) && !is.null(file)) {
    fit <- write_ideal_adaptor_stanfit(x = fit, file = file, compress = file_compress)
  }

  return(fit)
}
