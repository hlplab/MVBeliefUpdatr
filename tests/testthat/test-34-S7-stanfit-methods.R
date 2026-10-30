
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

  expect_equal(get_family_name(fit_1cue), "NIX")
  if (!is.null(fit_3cue)) {
    expect_equal(get_family_name(fit_3cue), "NIW")
  }
  expect_equal(get_family_name(get_staninput(fit_1cue)), "NIX")
})

test_that("summary.MVBU_Stanfit displays fixed parameters in @fixed", {
  fit_fixed <- fit_1cue
  fixed_spec <- c(
    fix_parameter(0.05, class = "lapse_rate"),
    fix_parameter(c(b = 0.5, p = 0.5), class = "category_prior")
  )
  fit_fixed@metadata$parameters <- list(
    fixed = list(
      original = fixed_spec,
      transformed = fixed_spec
    ),
    fitted = list("m", "S"),
    priors = list(
      original = list(),
      transformed = list()
    )
  )
  sum_fixed <- suppressWarnings(summary(fit_fixed))
  expect_true(S7::S7_inherits(sum_fixed, Summary_MVBU_Stanfit))
  expect_true(!is.null(sum_fixed@fixed))
  expect_equal(nrow(sum_fixed@fixed), 2L)
  expect_true("lapse_rate" %in% sum_fixed@fixed$class)
  expect_true("category_prior" %in% sum_fixed@fixed$class)
  # Fixed parameters should not be in fitted table
  expect_false("lapse_rate" %in% sum_fixed@fitted$Parameter)
})

test_that(
  "get_transform_information works for MVBU_Stanfit and MVBU_StanfitPosterior",
  {
    ti_fit <- get_transform_information(fit_1cue)
    expect_true(S7::S7_inherits(ti_fit, MVBU_TransformInformation))

    post_1cue <- as_MVBU_stanfit_posterior(fit_1cue)
    expect_true(S7::S7_inherits(post_1cue, MVBU_StanfitPosterior))
    ti_post <- get_transform_information(post_1cue)
    expect_identical(ti_post, ti_fit)

    custom_ti <- MVBU_TransformInformation()
    post_custom <- MVBU_StanfitPosterior(
      draws = list(),
      metadata = list(),
      data = list(
        original = tibble::tibble(),
        transformed = tibble::tibble()
      ),
      transform_information = custom_ti
    )
    expect_identical(get_transform_information(post_custom), custom_ti)

    expect_error(
      get_transform_information("invalid"),
      "x must be an MVBU_Stanfit, MVBU_StanfitPosterior"
    )
  }
)

test_that("get_data, get_parameters, and get_stan_priors space works", {
  data_orig <- get_data(fit_1cue, space = "original")
  data_trans <- get_data(fit_1cue, space = "transformed")
  expect_s3_class(data_orig, "data.frame")
  expect_s3_class(data_trans, "data.frame")

  exp_orig <- get_exposure_data(fit_1cue, space = "original")
  exp_trans <- get_exposure_data(fit_1cue, space = "transformed")
  expect_s3_class(exp_orig, "data.frame")
  expect_s3_class(exp_trans, "data.frame")

  test_orig <- get_test_data(fit_1cue, space = "original")
  test_trans <- get_test_data(fit_1cue, space = "transformed")
  expect_s3_class(test_orig, "data.frame")
  expect_s3_class(test_trans, "data.frame")

  pars_orig <- get_parameters(fit_1cue, space = "original")
  pars_trans <- get_parameters(fit_1cue, space = "transformed")
  expect_true(is.list(pars_orig))
  expect_true("fixed" %in% names(pars_orig))
  expect_true("fitted" %in% names(pars_orig))
  expect_false("prior" %in% names(pars_orig))
  expect_false("priors" %in% names(pars_orig))

  expect_true("fixed" %in% names(pars_trans))
  expect_true("fitted" %in% names(pars_trans))
  expect_false("prior" %in% names(pars_trans))
  expect_false("priors" %in% names(pars_trans))

  priors_orig <- get_stan_priors(fit_1cue, space = "original")
  priors_trans <- get_stan_priors(fit_1cue, space = "transformed")
  expect_true(!is.null(priors_orig))
  expect_true(!is.null(priors_trans))
})



