#' @include S7-stanfit-utils.R
#' @include S7-prior.R
NULL

#' Fit ideal adaptor
#'
#' Infers a prior and posterior distribution of ideal adaptors from the input
#' data using Stan. Currently, three types of ideal adaptor models are
#' available, each using a conjugate prior over Gaussian or multivariate
#' Gaussian categories:
#'
#' \describe{
#'   \item{`NIW_ideal_adaptor`}{A Normal-Inverse-Wishart (NIW) prior over the
#'     ideal adaptor. This is the default model with multivariate Gaussian
#'     categories. Accepts univariate and multivariate input, though the NIX
#'     model should be faster for univariate input.}
#'   \item{`NIX_ideal_adaptor`}{A Normal-Inverse-Chisquare (NIX) prior over
#'     the ideal adaptor with univariate Gaussian categories. Accepts only
#'     univariate input.}
#'   \item{`MNIX_ideal_adaptor`}{Separate NIXs for each of multiple cues that
#'     are integrated over during categorization assuming ideal cue weights
#'     based on the relative informativity of each cue.}
#' }
#'
#' @param stanfit_input An object of class
#'   \code{\link{IdealAdaptorStanfitInput}} containing the prepared data and
#'   transform metadata.
#' @param stanmodel Character string naming the Stan model to use. Defaults to
#'   \code{"NIW_ideal_adaptor"}.
#' @param fixed_parameters An optional \code{\link{MVBU_FixedParameter}} object
#'   specifying fixed cognitive parameters created via
#'   \code{\link{fix_parameter}}. Default: \code{NULL}.
#' @param priors An optional \code{\link{MVBU_Prior}} object specifying prior
#'   distributions created via \code{\link{prior}}. Default: \code{NULL}.
#' @param backend Character string naming the package to use as the backend for
#'   fitting the Stan model. Options are \code{"rstan"} (the default) or
#'   \code{"cmdstanr"}. Details on the \pkg{rstan} and \pkg{cmdstanr} packages
#'   are available at \url{https://mc-stan.org/rstan/} and
#'   \url{https://mc-stan.org/cmdstanr/}, respectively.
#' @param file Either NULL or a character string. In the latter case, the fitted
#'   model object is saved via saveRDS in a file named after the string supplied
#'   in file. The .rds extension is added automatically. If the file already
#'   exists, the model from that file will be loaded and returned instead of
#'   refitting the model. As existing files won't be overwritten, you have to
#'   manually remove the file in order to refit and save the model under an
#'   existing file name. The file name is stored in the
#'   \code{ideal_adaptor_stanfit} object for later use. (default: \code{NULL})
#' @param file_compress Logical or a character string, specifying one of the
#'   compression algorithms supported by \code{\link{saveRDS}}. If the
#'   \code{file} argument is provided, this compression will be used when saving
#'   the fitted model object.
#' @param file_refit Modifies when the fit stored via the \code{file} argument
#'   is re-used. For \code{"never"} (default) the fit is always loaded if it
#'   exists and fitting is skipped. For \code{"always"} the model is always
#'   refitted. If set to \code{"on_change"}, model will be refit if data passed
#'   to Stan differ from what is stored in the file. Refit will not be triggered
#'   for changes in additional parameters of the fit (e.g., initial values,
#'   number of iterations, control arguments, ...).
#' @param silent Verbosity level between \code{0} and \code{2}. If \code{1} (the
#'   default), most of the informational messages of compiler and sampler are
#'   suppressed. If \code{2}, even more messages are suppressed. The actual
#'   sampling progress is still printed. Set \code{refresh = 0} to turn this off
#'   as well. If using \code{backend = "rstan"} you can also set
#'   \code{open_progress = FALSE} to prevent opening additional progress bars.
#' @param threads Optional number of threads per chain (or list containing
#'   \code{threads}) for multithreading within Stan. Default: \code{NULL}.
#' @param stan_model_args A \code{list} of further arguments passed to
#'   [rstan::stan_model] for \code{backend = "rstan"} or to
#'   \code{cmdstanr::cmdstan_model} for \code{backend = "cmdstanr"}, which
#'   allows to change how models are compiled.
#' @param rename For internal use only.
#' @param ... Additional parameters are passed to [rstan::sampling], including
#'   the following:
#' @param stanvars Optional \code{stanvars} object. Default: \code{NULL}.
#' @param save_pars Storage for parameter names to save. Default: \code{NULL}.
#' @param chains Number of Markov chains. Default: 4.
#' @param iter Total number of iterations per chain. Default: 2000.
#' @param warmup Number of warmup (burnin) iterations per chain. Default: 1000.
#' @param init Initial values specification. Default: \code{"random"}.
#' @param control Named list of control parameters passed to Stan sampler.
#'   Default: \code{NULL}.
#'
#' @return An object of class \code{ideal_adaptor_stanfit} with the fitted
#'   stan model.
#'
#' @details
#'   In interpreting the inferred parameters, it should be kept in mind that
#'   the \emph{inferred} scatter matrix (e.g., S_0 for the NIW model) includes
#'   variability from internal perceptual and/or external environmental noise,
#'   \emph{in addition} to the motor noise that is reflected in production
#'   data. This also implies that, \strong{if \code{Sigma_0} is provided
#'   by the user it should be arguably convolved with an estimate of
#'   perceptual noise}. This is particularly important if the data you're
#'   fitting contains test phases without exposure (e.g., pre-exposure tests).
#'
#'   Use \code{methods(class = "ideal_adaptor_stanfit")} for an overview on
#'   available methods.
#'
#' @seealso \code{\link{IdealAdaptorStanfit}} for information about ideal
#'   adaptor stanfit objects, \code{\link{get_draws}} to draw samples from the
#'   stanfit.
#'
#' @importFrom rstan nlist
#' @export
fit_ideal_adaptor <- function(
  stanfit_input,
  stanmodel = "NIW_ideal_adaptor",
  fixed_parameters = NULL,
  priors = NULL,
  file = NULL,
  file_refit = "never",
  file_compress = TRUE,
  # included for later use
  stanvars = NULL,
  backend = "rstan",
  save_pars = NULL,
  chains = 4,
  iter = 2000,
  warmup = 1000,
  init = "random",
  control = NULL,
  silent = 1,
  threads = NULL,
  stan_model_args = list(),
  # Stuff to be deprecated in the future
  rename = TRUE,
  ...
) {
  backend <- match.arg(.as_one_character(backend), c("rstan", "cmdstanr"))
  file_refit <- match.arg(file_refit, .file_refit_options())

  # optionally load ideal_adaptor_stanfit from file
  if (!is.null(file) && file_refit == "never") {
    fit <- read_stanfit(file)
    if (!is.null(fit)) {
      if (silent == 0) message("Loading existing model from file.")
      return(fit)
    }
  }

  .assert_true(
    S7::S7_inherits(stanfit_input, IdealAdaptorStanfitInput),
    msg = "stanfit_input must be an IdealAdaptorStanfitInput object."
  )

  if (is.null(stanmodel)) {
    stanmodel <- "NIW_ideal_adaptor"
  }

  .assert_that(
    is.character(stanmodel) && length(stanmodel) == 1L && !is.na(stanmodel) &&
      !is.null(stanmodels[[stanmodel]]),
    msg = paste(
      "The specified stanmodel does not exist. Allowable models include:",
      paste(names(stanmodels), collapse = ", ")
    )
  )

  data <- list(
    original = get_data(stanfit_input, space = "original"),
    transformed = get_data(stanfit_input, space = "transformed")
  )
  transform_information <- stanfit_input@transform_information
  cues <- get_cue_labels(stanfit_input)
  categories <- get_category_labels(stanfit_input)
  groups <- get_group_labels(stanfit_input, include_prior = FALSE)

  spec <- .validate_and_transform_model_spec(
    stanmodel = stanmodel,
    cues = cues,
    categories = categories,
    groups = groups,
    fixed_parameters = fixed_parameters,
    priors = priors,
    transform = transform_information,
    control = if (is.list(control)) control else list()
  )

  staninput <- .new_staninput(
    stanfit_input = stanfit_input,
    stanmodel = stanmodel,
    fixed_parameters_transformed = spec$fixed_transformed,
    priors_transformed = spec$priors_transformed,
    control = if (is.list(control)) control else list()
  )

  # Check whether model actually needs to be refit
  if (!is.null(file) && file_refit == "on_change") {
    x_from_file <- read_stanfit(file)
    if (!is.null(x_from_file)) {
      needs_refit <- .stanfit_needs_refit(
        x_from_file,
        current_version = get_current_versions(),
        data = data$original,
        staninput = staninput,
        silent = (silent >= 1)
      )
      if (!needs_refit) {
        if (silent == 0) {
          message("No refitting needed. Loading existing model from file.")
        }
        return(x_from_file)
      }
    }
  }

  fam <- .to_family(stanmodel)

  metadata <- if (
    !is.null(stanfit_input@metadata) && length(stanfit_input@metadata) > 0
  ) {
    stanfit_input@metadata
  } else {
    list(label_information = get_labels(stanfit_input))
  }
  metadata$model_family <- fam
  metadata$parameters <- list(
    fixed = list(
      original = spec$fixed_untransformed,
      transformed = spec$fixed_transformed
    ),
    fitted = as.list(
      setdiff(
        c("mu_0", "Sigma_0", "kappa_0", "nu_0"),
        names(spec$fixed_untransformed)
      )
    ),
    priors = list(
      original = spec$priors_untransformed,
      transformed = spec$priors_transformed$table
    )
  )

  constructor <- .get_ideal_adaptor_stanfit_constructor(
    staninput = staninput,
    stanmodel = stanmodel
  )
  fit <- constructor(
    data = data,
    staninput = staninput,
    stanvars = stanvars,
    save_pars = save_pars,
    backend = backend,
    stan_args = .nlist(init, silent, control, threads, stan_model_args),
    transform_information = transform_information,
    file = file,
    metadata = metadata
  )

  # Check that staninput has at least two categories (fitting with one
  # category makes no sense)
  .assert_true(
    get_staninput(fit)@values$M >= 2,
    msg = "staninput must have at least two categories."
  )

  stanfit <- NULL
  if (chains > 0 & iter > 0) {
    # Parameters *not* to store (used with include = F below)
    # Future file reduction could be achieved via the shredder package for
    # stanfit post-processing (https://github.com/yonicd/shredder)
    exclude_pars <- c(
      "lapse_rate_param",
      "m_0_param", "m_0_tau", "m_0_tau_param",
      "m_0_L_omega", "m_0_L_omega_param",
      "tau_0_param", "L_omega_0_param", "L_S_0", "L_S_n", "L_t_scale",
      "p_test_conj", "log_p_test_conj",
      "z_mean_exposure", "z_css_exposure", "z_test",
      "mu_0_1d_data", "Sigma_0_1d_data"
    )

    stanfit <- .fit_model(
      stanmodel = stanmodel,
      backend = backend,
      data = get_staninput(fit)@values,
      pars = exclude_pars,
      include = FALSE,
      chains = chains,
      iter = iter,
      warmup = warmup,
      init = init,
      control = control,
      threads = threads,
      silent = silent,
      stan_model_args = stan_model_args,
      ...
    )

    .assert_contains_draws(stanfit)
    if (requireNamespace("tidybayes", quietly = TRUE)) {
      category_labels <- get_category_labels(stanfit_input)
      group_labels <- get_group_labels(stanfit_input)
      cue_labels <- get_cue_labels(stanfit_input)
      stanfit <- tidybayes::recover_types(
        stanfit,
        list(
          category = factor(category_labels, levels = category_labels),
          group = factor(group_labels, levels = group_labels),
          cue = factor(cue_labels, levels = cue_labels),
          cue2 = factor(cue_labels, levels = cue_labels)
        )
      )
    }
    fit <- set_stanfit(fit, stanfit)
    sampled_pars <- .extract_fitted_parameter_names_from_stanfit(stanfit)
    fixed_names <- names(get_parameters(fit, space = "original")$fixed)
    fit@metadata$parameters$fitted <-
      as.list(setdiff(sampled_pars, fixed_names))
  } else if (!silent) {
    message("No sampling requested. Returning empty model object.")
  }

  if (!is.null(fit) && !is.null(file)) {
    fit <- write_stanfit(x = fit, file = file, compress = file_compress)
  }

  return(fit)
}

# -------------------------------------------------------------------------
# Backend Model Fitting Dispatchers
# -------------------------------------------------------------------------

#' Fit a Stan model using the requested backend
#'
#' @param stanmodel Character string naming the Stan model.
#' @param backend Character string naming the backend: \code{"rstan"} or
#'   \code{"cmdstanr"}.
#' @param data Named list of Stan input data.
#' @param pars Character vector of parameter names to include/exclude.
#' @param include Logical; if \code{FALSE}, parameters in \code{pars} are
#'   excluded.
#' @param chains Integer number of chains.
#' @param iter Total iterations per chain.
#' @param warmup Warmup iterations per chain.
#' @param init Initial values specification.
#' @param control List of sampler control options.
#' @param threads Thread specification (count or list).
#' @param silent Verbosity level (0 to 2).
#' @param stan_model_args List of compilation options.
#' @param ... Additional arguments passed to the backend sampler.
#' @return An object of class \code{\link[rstan:stanfit-class]{stanfit}}.
#' @keywords internal
#' @noRd
.fit_model <- function(
  stanmodel,
  backend = "rstan",
  data,
  pars = character(0),
  include = FALSE,
  chains = 4,
  iter = 2000,
  warmup = 1000,
  init = "random",
  control = NULL,
  threads = NULL,
  silent = 1,
  stan_model_args = list(),
  ...
) {
  fit_fun <- get(paste0(".fit_model_", backend), mode = "function")
  fit_fun(
    stanmodel = stanmodel,
    data = data,
    pars = pars,
    include = include,
    chains = chains,
    iter = iter,
    warmup = warmup,
    init = init,
    control = control,
    threads = threads,
    silent = silent,
    stan_model_args = stan_model_args,
    ...
  )
}

#' Fit a Stan model using rstan
#'
#' @inheritParams .fit_model
#' @keywords internal
#' @noRd
.fit_model_rstan <- function(
  stanmodel,
  data,
  pars = character(0),
  include = FALSE,
  chains = 4,
  iter = 2000,
  warmup = 1000,
  init = "random",
  control = NULL,
  threads = NULL,
  silent = 1,
  stan_model_args = list(),
  ...
) {
  threads_spec <- .validate_threads(threads)
  use_threading <- isTRUE(threads_spec$threads > 1L)
  use_precompiled <- (stanmodel %in% names(stanmodels)) &&
    length(stan_model_args) == 0L &&
    !use_threading

  if (use_precompiled) {
    model_obj <- stanmodels[[stanmodel]]
  } else {
    stan_file <- .get_stan_file(stanmodel)
    compile_args <- c(
      list(
        model = stan_file, backend = "rstan",
        threads = threads_spec, silent = silent
      ),
      stan_model_args
    )
    model_obj <- .do_call(.compile_model, compile_args)
  }

  if (use_threading) {
    if (utils::packageVersion("rstan") >= "2.26") {
      threads_per_chain_def <- rstan::rstan_options("threads_per_chain")
      on.exit(
        rstan::rstan_options(threads_per_chain = threads_per_chain_def),
        add = TRUE
      )
      rstan::rstan_options(threads_per_chain = threads_spec$threads)
    } else {
      .stop(
        "Threading is not supported by backend 'rstan' version ",
        utils::packageVersion("rstan"), "."
      )
    }
  }

  sampling_args <- list(
    object = model_obj,
    data = data,
    check_data = TRUE,
    pars = pars,
    include = include,
    chains = chains,
    iter = iter,
    warmup = warmup,
    init = init,
    control = control,
    show_messages = (silent < 1)
  )
  dots <- list(...)
  sampling_args[names(dots)] <- dots

  .do_call(rstan::sampling, sampling_args)
}

#' Fit a Stan model using cmdstanr
#'
#' @inheritParams .fit_model
#' @keywords internal
#' @noRd
.fit_model_cmdstanr <- function(
  stanmodel,
  data,
  pars = character(0),
  include = FALSE,
  chains = 4,
  iter = 2000,
  warmup = 1000,
  init = "random",
  control = NULL,
  threads = NULL,
  silent = 1,
  stan_model_args = list(),
  ...
) {
  if (!requireNamespace("cmdstanr", quietly = TRUE)) {
    .stop("Package 'cmdstanr' is required to use backend 'cmdstanr'.")
  }

  stan_file <- .get_stan_file(stanmodel)
  compile_args <- c(
    list(
      model = stan_file, backend = "cmdstanr",
      threads = threads, silent = silent
    ),
    stan_model_args
  )
  model_obj <- .do_call(.compile_model, compile_args)

  sample_args <- list(
    data = data,
    chains = as.integer(chains),
    iter_sampling = as.integer(iter - warmup),
    iter_warmup = as.integer(warmup),
    show_messages = (silent < 2),
    show_exceptions = (silent == 0)
  )

  if (identical(init, "random")) {
    sample_args$init <- NULL
  } else if (identical(init, "0") || identical(init, 0)) {
    sample_args$init <- 0
  } else {
    sample_args$init <- init
  }

  if (is.list(control)) {
    sample_args[names(control)] <- control
  }

  threads_spec <- .validate_threads(threads)
  if (isTRUE(threads_spec$threads > 1L)) {
    sample_args$threads_per_chain <- threads_spec$threads
  }

  dots <- list(...)
  sample_args[names(dots)] <- dots

  cs_fit <- .do_call(model_obj$sample, sample_args)

  exclude_pars <- if (!include) pars else character(0)
  .read_cmdstan_csv_as_stanfit(
    files = cs_fit$output_files(),
    exclude_pars = exclude_pars,
    model_name = stanmodel
  )
}
