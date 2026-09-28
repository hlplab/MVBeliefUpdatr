# ==============================================================================
# Global Configuration for Stanfit Test Fixtures
# ==============================================================================
# `refit_stanfit_models`: Controls whether full Stanfit test fixtures are
# refitted during test execution.
#
# Default: FALSE.
# - When FALSE: Tests load existing pre-fitted fixtures from disk with
#   `file_refit = "never"`. Tests will NOT be skipped if a fixture is missing
#   or if an outdated cached model causes downstream assertion failures;
#   instead, tests will fail explicitly, signalling to the developer that the
#   fixtures must be refitted.
# - When TRUE: Full models are refitted using `file_refit = "on_change"`.
#   Fitting full models (1-cue NIX, 2-cue MNIX, 3-cue NIW) can take substantial
#   time (an hour or longer the first time Stan code or dependencies
#   change). Once refitted and saved to disk, set `refit_stanfit_models <- FALSE`
#   (or unset the environment variable) so that subsequent test runs execute
#   rapidly by loading the cached fixtures.
#
# Developers can toggle this variable below or set the environment variable:
#   Sys.setenv(MVBU_REFIT_STANFIT = "true")
# ==============================================================================
refit_stanfit_models <- isTRUE(
  as.logical(Sys.getenv("MVBU_REFIT_STANFIT", "FALSE"))
) || FALSE

.get_ideal_adaptor_fit_model_dir <- function() {
  testthat::test_path("models")
}

example_stanfit_path <- function(
  example,
  stanmodel = "NIW_ideal_adaptor",
  seed = 42,
  control = control_staninput()
) {
  transform_type <- control$transform_type
  file.path(
    .get_ideal_adaptor_fit_model_dir(),
    paste0(
      "example-stanfit-",
      paste(
        c(
          if (!is.null(stanmodel)) stanmodel else "",
          example,
          if (!is.null(transform_type)) transform_type else "",
          seed
        ),
        collapse = "-"
      ),
      ".rds"
    )
  )
}

get_example_stanfit <- function(
  example = 1L,
  silent = 2,
  refresh = 0,
  seed = 42L,
  verbose = FALSE,
  file_refit = if (isTRUE(refit_stanfit_models)) "on_change" else "never",
  stanmodel = "NIW_ideal_adaptor",
  lapse_rate = NULL,
  mu_0 = NULL,
  Sigma_0 = NULL,
  control = control_staninput(),
  transform_type = NULL,
  filename = NULL,
  file = NULL,
  ...
) {
  if (!is.null(transform_type)) {
    control$transform_type <- transform_type
  }
  if (is.null(filename)) {
    filename <- if (!is.null(file)) {
      file
    } else {
      example_stanfit_path(
        example,
        stanmodel = stanmodel,
        seed = seed,
        control = control
      )
    }
  }

  model_family <- sub("_ideal_adaptor$", "", stanmodel)
  example_ideal_adaptor_stanfit(
    model_family = model_family,
    n_cues = example,
    seed = seed,
    staninput_control = control,
    file = filename,
    file_refit = file_refit,
    silent = silent,
    refresh = refresh,
    ...
  )
}

get_full_stanfit <- function(
  family = c("NIX", "MNIX", "NIW"),
  file_refit = if (isTRUE(refit_stanfit_models)) "on_change" else "never",
  ...
) {
  family <- match.arg(family)
  switch(
    family,
    NIX = get_example_stanfit(
      example = 1L,
      stanmodel = "NIX_ideal_adaptor",
      file_refit = file_refit,
      ...
    ),
    MNIX = get_example_stanfit(
      example = 2L,
      stanmodel = "MNIX_ideal_adaptor",
      file_refit = file_refit,
      ...
    ),
    NIW = get_example_stanfit(
      example = 3L,
      stanmodel = "NIW_ideal_adaptor",
      file_refit = file_refit,
      ...
    )
  )
}

get_minimal_stanfit <- function(
  family = c("NIX", "MNIX", "NIW"),
  n_cues = NULL,
  file_refit = "on_change",
  ...
) {
  family <- match.arg(family)
  if (is.null(n_cues)) {
    n_cues <- if (family == "NIX") 1L else 2L
  }
  filename <- switch(
    family,
    NIX = "minimal-nix_ideal_adaptor-NIX.rds",
    MNIX = "minimal-mnix_ideal_adaptor-MNIX-2cue.rds",
    NIW = if (n_cues == 1L) {
      "minimal-niw_ideal_adaptor-NIW-1cue.rds"
    } else {
      "minimal-niw_ideal_adaptor-NIW-2cue.rds"
    }
  )
  model_path <- testthat::test_path("models", filename)
  example_ideal_adaptor_stanfit(
    model_family = family,
    n_cues = n_cues,
    seed = 42L,
    file = model_path,
    file_refit = file_refit,
    chains = 1,
    iter = 500,
    warmup = 250,
    refresh = 0,
    control = list(adapt_delta = 0.99, max_treedepth = 12),
    ...
  )
}

expect_staninput_structure <- function(
  input,
  expected_class,
  required_names,
  forbidden_names = character()
) {
  expect_true(S7::S7_inherits(input, IdealAdaptorStanfitInput))
  expect_true(S7::S7_inherits(input@staninput, expected_class))
  expect_true(S7::S7_inherits(input@staninput, IdealAdaptorStaninput))
  expect_true(is.list(input@staninput@values))
  expect_true(all(required_names %in% names(input@staninput@values)))
  expect_false(any(forbidden_names %in% names(input@staninput@values)))
}

run_fixed_param_stan_program <- function(stan_file, input, model_name) {
  skip_if_not_installed("rstan")

  model_key <- sub("_compat$", "", model_name)
  model <- NULL
  if (
    exists(
      "stanmodels",
      envir = asNamespace("MVBeliefUpdatr"),
      inherits = FALSE
    )
  ) {
    pkg_models <- get("stanmodels", envir = asNamespace("MVBeliefUpdatr"))
    if (!is.null(pkg_models[[model_name]])) {
      model <- pkg_models[[model_name]]
    } else if (!is.null(pkg_models[[model_key]])) {
      model <- pkg_models[[model_key]]
    } else if (!is.null(pkg_models[[paste0(model_key, "_ideal_adaptor")]])) {
      model <- pkg_models[[paste0(model_key, "_ideal_adaptor")]]
    }
  }

  if (is.null(model)) {
    model <- rstan::stan_model(
      file = stan_file,
      model_name = model_name,
      auto_write = FALSE
    )
  }

  fit <- rstan::sampling(
    object = model,
    data = input@staninput@values,
    chains = 1,
    iter = 1,
    warmup = 0,
    refresh = 0,
    seed = 123,
    algorithm = "Fixed_param"
  )

  expect_s4_class(fit, "stanfit")
}
