test_that("fit_ideal_adaptor fits minimal shifted-prior examples for the supported models", {
  skip_if_not_installed("rstan")

  for (spec in list(
    list(name = "NIX", cues = "VOT", stanmodel = "NIX_ideal_adaptor"),
    list(name = "NIW_1cue", cues = "VOT", stanmodel = "NIW_ideal_adaptor"),
    list(name = "NIW_2cue", cues = c("VOT", "f0_semitones"), stanmodel = "NIW_ideal_adaptor"),
    list(name = "MNIX", cues = c("VOT", "f0_semitones"), stanmodel = "MNIX_ideal_adaptor")
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
      iter = 450,
      warmup = 225,
      refresh = 0,
      control = list(adapt_delta = 0.95),
      stanmodel = spec$stanmodel
    )

    expect_true(S7::S7_inherits(fit, IdealAdaptorStanfit), info = paste("fit object is not an ideal adaptor fit for", spec$name))
    expect_s4_class(get_stanfit(fit), "stanfit")
    expect_true(.contains_draws(get_stanfit(fit)), info = paste("fit does not contain posterior draws for", spec$name))
  }
})
