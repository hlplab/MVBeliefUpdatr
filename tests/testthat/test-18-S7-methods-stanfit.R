
skip_if_not(file.exists(example_stanfit_path(1, stanmodel = "NIW_ideal_adaptor")), "cached example stanfit not generated")
fit_1cue <- get_example_stanfit(1, stanmodel = "NIW_ideal_adaptor", file_refit = "never")

test_that("Test for single cue", {
  expect_no_error(suppressWarnings(summary(fit_1cue)))
  expect_no_error(suppressWarnings(summary(fit_1cue, only_prior = TRUE)))
  expect_no_error(suppressWarnings(summary(fit_1cue, include_transformed_pars = TRUE)))
  expect_true(S7::S7_inherits(suppressWarnings(summary(fit_1cue)), Summary_MVBU_Stanfit) || is.data.frame(suppressWarnings(summary(fit_1cue))))
})

fit_3cue <- tryCatch(
  get_example_stanfit(3, stanmodel = "NIW_ideal_adaptor", file_refit = "never"),
  error = function(e) NULL
)

test_that("Test for multiple cues", {
  skip_if(is.null(fit_3cue), "3-cue model not yet available or needs refit")
  expect_no_error(suppressWarnings(summary(fit_3cue)))
  expect_no_error(suppressWarnings(summary(fit_3cue, only_prior = TRUE)))
  expect_no_error(suppressWarnings(summary(fit_3cue, include_transformed_pars = TRUE)))
  expect_true(S7::S7_inherits(suppressWarnings(summary(fit_3cue)), Summary_MVBU_Stanfit) || is.data.frame(suppressWarnings(summary(fit_3cue))))
})

test_that("Stanfit diagnostic and posterior methods", {
  rh <- rhat(fit_1cue)
  expect_true(is.numeric(rh))
  expect_true(length(rh) > 0)
  expect_true(all(rh >= 0, na.rm = TRUE))

  neff <- neff_ratio(fit_1cue)
  expect_true(is.numeric(neff))
  expect_true(length(neff) > 0)

  lp <- log_posterior(fit_1cue)
  expect_true(is.data.frame(lp))
  expect_true("Value" %in% names(lp))

  np <- nuts_params(fit_1cue)
  expect_true(is.data.frame(np))
  expect_true("Parameter" %in% names(np))

  cp <- control_params(fit_1cue)
  expect_true(is.list(cp))
  expect_true("adapt_delta" %in% names(cp) || "max_treedepth" %in% names(cp))

  expect_true(inherits(posterior::as_draws(fit_1cue), "draws"))
  expect_true(inherits(posterior::as_draws_df(fit_1cue), "draws_df"))
  expect_true(inherits(posterior::as_draws_array(fit_1cue), "draws_array"))
  expect_true(inherits(posterior::as_draws_matrix(fit_1cue), "draws_matrix"))
  expect_true(inherits(posterior::as_draws_list(fit_1cue), "draws_list"))
  expect_true(inherits(posterior::as_draws_rvars(fit_1cue), "draws_rvars"))

  expect_equal(get_model_type(fit_1cue), "NIW_ideal_adaptor")
  if (!is.null(fit_3cue)) {
    expect_equal(get_model_type(fit_3cue), "NIW_ideal_adaptor")
  }
  expect_equal(get_model_type(get_staninput(fit_1cue)), "NIW_ideal_adaptor")
})

