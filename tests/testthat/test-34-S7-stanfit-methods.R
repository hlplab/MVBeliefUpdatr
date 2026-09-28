
fit_1cue <- get_full_stanfit("NIX")

test_that("Test for single cue", {
  sum_1cue <- suppressWarnings(summary(fit_1cue))
  expect_no_error(suppressWarnings(summary(fit_1cue)))
  expect_no_error(suppressWarnings(summary(fit_1cue, only_prior = TRUE)))
  expect_no_error(
    suppressWarnings(summary(fit_1cue, include_transformed_pars = TRUE))
  )
  expect_true(S7::S7_inherits(sum_1cue, Summary_MVBU_Stanfit))
  expect_true("Cue" %in% names(sum_1cue@fitted))
  expect_false("Cue1" %in% names(sum_1cue@fitted))
  expect_false("Cue2" %in% names(sum_1cue@fitted))
  m_rows <- sum_1cue@fitted[sum_1cue@fitted$Parameter == "m", ]
  expect_true(all(!is.na(m_rows$Cue) & m_rows$Cue != ""))
  S_rows <- sum_1cue@fitted[sum_1cue@fitted$Parameter == "S", ]
  expect_true(all(!is.na(S_rows$Cue) & S_rows$Cue != ""))
})

fit_3cue <- get_full_stanfit("NIW")

test_that("Test for multiple cues", {
  sum_3cue <- suppressWarnings(summary(fit_3cue))
  expect_no_error(suppressWarnings(summary(fit_3cue)))
  expect_no_error(suppressWarnings(summary(fit_3cue, only_prior = TRUE)))
  expect_no_error(
    suppressWarnings(summary(fit_3cue, include_transformed_pars = TRUE))
  )
  expect_true(S7::S7_inherits(sum_3cue, Summary_MVBU_Stanfit))
  expect_true("Cue1" %in% names(sum_3cue@fitted))
  expect_true("Cue2" %in% names(sum_3cue@fitted))
  expect_false("Cue" %in% names(sum_3cue@fitted))
})

test_that("MNIX summary has Cue and no Cue2", {
  fit_mnix <- get_full_stanfit("MNIX")
  sum_mnix <- suppressWarnings(summary(fit_mnix))
  expect_true(S7::S7_inherits(sum_mnix, Summary_MVBU_Stanfit))
  expect_true("Cue" %in% names(sum_mnix@fitted))
  expect_false("Cue1" %in% names(sum_mnix@fitted))
  expect_false("Cue2" %in% names(sum_mnix@fitted))
  if (!is.null(sum_mnix@high_rhats) && nrow(sum_mnix@high_rhats) > 0) {
    out_mnix <- capture.output(print(sum_mnix))
    expect_true(any(grepl("Parameters with Rhats > 1.05:", out_mnix)))
  }
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

  expect_equal(get_model_type(fit_1cue), "NIX_ideal_adaptor")
  if (!is.null(fit_3cue)) {
    expect_equal(get_model_type(fit_3cue), "NIW_ideal_adaptor")
  }
  expect_equal(get_model_type(get_staninput(fit_1cue)), "NIX_ideal_adaptor")
})

