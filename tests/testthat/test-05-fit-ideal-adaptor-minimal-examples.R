test_that("fit_ideal_adaptor fits minimal shifted-prior examples for the current NIX/NIW models", {
  skip_if_not_installed("rstan")

  for (spec in list(
    list(name = "NIX", cues = "VOT", stanmodel = "NIX_ideal_adaptor"),
    list(name = "NIW_1cue", cues = "VOT", stanmodel = "NIW_ideal_adaptor"),
    list(name = "NIW_2cue", cues = c("VOT", "f0_semitones"), stanmodel = "NIW_ideal_adaptor")
    # list(name = "MNIX", cues = c("VOT", "f0_semitones"), stanmodel = "MNIX_ideal_adaptor")
  )) {
    data <- build_shifted_prior_fit_example(
      cues = spec$cues,
      n_exposure_per_condition = 40L,
      n_test_per_condition = 30L,
      seed = 123L + match(spec$name, c("NIX", "NIW_1cue", "NIW_2cue", "MNIX"))
    )

    input <- new_ideal_adaptor_stanfit_input(
      exposure = data$exposure,
      test = data$test,
      cues = spec$cues,
      category = "category",
      response = "response",
      group = "group",
      group.unique = "Condition",
      control = control_staninput(transform_type = "identity"),
      stanmodel = spec$stanmodel
    )

    fit <- fit_ideal_adaptor(
      stanfit_input = input,
      chains = 1,
      iter = 600,
      warmup = 300,
      refresh = 0,
      control = list(adapt_delta = 0.99, max_treedepth = 12),
      stanmodel = spec$stanmodel
    )

    model_name <- paste0("minimal-", tolower(spec$stanmodel), "-", gsub("[^A-Za-z0-9]+", "-", spec$name))
    model_path <- save_ideal_adaptor_fit_model(fit, model_name)
    expect_true(file.exists(model_path), info = paste("model file was not written for", spec$name))

    reloaded_fit <- load_ideal_adaptor_fit_model(model_name)
    expect_true(S7::S7_inherits(reloaded_fit, IdealAdaptorStanfit), info = paste("reloaded model is not an ideal adaptor fit for", spec$name))

    expect_true(S7::S7_inherits(fit, IdealAdaptorStanfit), info = paste("fit object is not an ideal adaptor fit for", spec$name))
    expect_s4_class(get_stanfit(fit), "stanfit")
    expect_true(.contains_draws(get_stanfit(fit)), info = paste("fit does not contain posterior draws for", spec$name))
  }
})


test_that("fit_ideal_adaptor loads an existing model from file when file_refit is never", {
  skip_if_not_installed("rstan")

  model_name <- "minimal-nix_ideal_adaptor-nix"
  fit <- load_ideal_adaptor_fit_model(model_name)
  tmp_file <- tempfile(fileext = ".rds")
  saveRDS(fit, tmp_file)

  data <- build_shifted_prior_fit_example(
    cues = "VOT",
    n_exposure_per_condition = 20L,
    n_test_per_condition = 10L,
    seed = 321L
  )

  input <- new_ideal_adaptor_stanfit_input(
    exposure = data$exposure,
    test = data$test,
    cues = "VOT",
    category = "category",
    response = "response",
    group = "group",
    group.unique = "Condition",
    control = control_staninput(transform_type = "identity"),
    stanmodel = "NIX_ideal_adaptor"
  )

  reloaded_fit <- fit_ideal_adaptor(
    stanfit_input = input,
    file = tmp_file,
    file_refit = "never",
    chains = 0,
    iter = 0,
    warmup = 0,
    refresh = 0,
    stanmodel = "NIX_ideal_adaptor"
  )

  expect_true(S7::S7_inherits(reloaded_fit, IdealAdaptorStanfit))
  # stanfit@.MISC holds rstan's C++ module pointer, which is rebuilt on every
  # deserialization, so compare the recovered content rather than the whole object.
  expect_equal(reloaded_fit@stanfit@sim$samples, fit@stanfit@sim$samples)
  expect_identical(reloaded_fit@stanfit@model_name, fit@stanfit@model_name)
  expect_identical(reloaded_fit@stanfit@model_pars, fit@stanfit@model_pars)
})

test_that("recover_types works on stanfit objects nested in S7 IdealAdaptorStanfit objects", {
  skip_if_not_installed("rstan")
  skip_if_not_installed("tidybayes")

  model_name <- "minimal-nix_ideal_adaptor-nix"
  fit <- load_ideal_adaptor_fit_model(model_name)

  recovered_stanfit <- tidybayes::recover_types(get_stanfit(fit))
  recovered_fit <- set_stanfit(fit, recovered_stanfit)

  expect_true(S7::S7_inherits(recovered_fit, IdealAdaptorStanfit))
  expect_s4_class(get_stanfit(recovered_fit), "stanfit")
  expect_true(!is.null(attr(get_stanfit(recovered_fit), "tidybayes_constructors")))
  expect_true(is.function(get_constructor(recovered_fit, "group")))
  expect_equal(get_staninput_variable_levels(recovered_fit, "group"), get_group_levels(recovered_fit))
})