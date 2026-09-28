fit <- get_full_stanfit("NIX")

test_that("add draws - input check (1 cue)", {
  expect_true(is_tibble(get_draws(fit, groups = "prior")))
  expect_true(is_tibble(get_draws(fit, groups = "left_shifted")))
  expect_true(is_tibble(get_draws(fit, groups = c("prior", "left_shifted"))))
  expect_true(is_tibble(get_draws(fit, groups = "prior", ndraws = 1)))
  expect_error(get_draws(fit, groups = "priors"))
  expect_error(get_draws(fit, groups = "prior", ndraws = c(1, 2)))
})

test_that("add draws - output check (1 cue)", {
  expect_equal(
    length(
      unique(get_draws(fit, groups = "prior", ndraws = 10, seed = 1)$.draw)
    ),
    10
  )
  expect_equal(
    nrow(get_draws(fit, groups = "prior", summarize = TRUE)),
    length(get_category_labels(fit))
  )
  expect_equal(
    names(get_draws(fit, groups = "prior", summarize = TRUE)),
    c(".chain", ".iteration", ".draw", "group", "category", "kappa", "nu", "m",
      "S", "lapse_rate", "Sigma_exp", "Sigma_marg")
  )
  expect_equal(
    names(get_draws(fit, groups = "prior", summarize = TRUE, nest = TRUE)),
    c(".chain", ".iteration", ".draw", "group", "category", "kappa", "nu", "m",
      "S", "lapse_rate", "Sigma_exp", "Sigma_marg")
  )
  expect_equal(
    names(get_draws(fit, groups = "prior", summarize = TRUE, nest = FALSE)),
    c(".chain", ".iteration", ".draw", "group", "category", "kappa", "nu",
      "cue", "cue2", "m", "S", "lapse_rate", "Sigma_exp", "Sigma_marg")
  )
})

fit <- get_full_stanfit("MNIX")

test_that("add draws - input check (2 cues)", {
  expect_true(is_tibble(get_draws(fit, groups = "prior")))
  expect_true(is_tibble(get_draws(fit, groups = "top_right")))
  expect_true(is_tibble(get_draws(fit, groups = c("prior", "top_right"))))
  expect_true(is_tibble(get_draws(fit, groups = "prior", ndraws = 1)))
  expect_error(get_draws(fit, groups = "priors"))
  expect_error(get_draws(fit, groups = "prior", ndraws = c(1, 2)))
})

test_that("add draws - output check (2 cues)", {
  expect_equal(
    length(
      unique(get_draws(fit, groups = "prior", ndraws = 10, seed = 1)$.draw)
    ),
    10
  )
  expect_equal(
    nrow(get_draws(fit, groups = "prior", summarize = TRUE)),
    length(get_category_labels(fit))
  )
  expect_equal(
    names(get_draws(fit, groups = "prior", summarize = TRUE)),
    c(".chain", ".iteration", ".draw", "group", "category", "kappa", "nu", "m",
      "S", "lapse_rate", "Sigma_exp", "Sigma_marg")
  )
  expect_equal(
    names(get_draws(fit, groups = "prior", summarize = TRUE, nest = TRUE)),
    c(".chain", ".iteration", ".draw", "group", "category", "kappa", "nu", "m",
      "S", "lapse_rate", "Sigma_exp", "Sigma_marg")
  )
  expect_equal(
    names(get_draws(fit, groups = "prior", summarize = TRUE, nest = FALSE)),
    c(".chain", ".iteration", ".draw", "group", "category", "kappa", "nu",
      "cue", "cue2", "m", "S", "lapse_rate", "Sigma_exp", "Sigma_marg")
  )
})

test_that("get exposure category statistic", {
  # Get error when *exposure* statistics is requested for *prior*
  expect_error(get_exposure_category_statistic(fit, groups = "prior"))
  expect_error(get_exposure_category_mean(fit, groups = "prior"))
  expect_error(get_exposure_category_cov(fit, groups = "prior"))
  # When prior is not requested
  # get_exposure_mean works for MNIX
  expect_true(is.vector(get_exposure_category_mean(fit, "/d/", "no_exposure")))
  expect_error(
    is.vector(get_exposure_category_mean(fit, "wrong", "no_exposure"))
  )
  expect_error(is.vector(get_exposure_category_mean(fit, "/d/", "wrong")))
  expect_true(
    is_tibble(
      get_exposure_category_mean(
        fit, c("/d/", "/t/"), c("no_exposure", "top_right")
      )
    )
  )
  # css, uss, cov for MNIX return matrices
  expect_true(is.matrix(get_exposure_category_css(fit, "/d/", "no_exposure")))
  expect_true(is.matrix(get_exposure_category_uss(fit, "/d/", "no_exposure")))
  expect_true(is.matrix(get_exposure_category_cov(fit, "/d/", "no_exposure")))
})

test_that("get expected category statistic", {
  expect_true(is.vector(get_expected_mu(fit, "/d/", "prior")))
  expect_error(is.vector(get_expected_mu(fit, "wrong", "prior")))
  expect_error(is.vector(get_expected_mu(fit, "/d/", "wrong")))
  expect_true(is.matrix(get_expected_sigma(fit, "/d/", "prior")))
  expect_error(is.matrix(get_expected_sigma(fit, "wrong", "prior")))
  expect_error(is.matrix(get_expected_sigma(fit, "/d/", "wrong")))
  expect_true(
    is_tibble(
      get_expected_sigma(fit, c("/d/", "/t/"), c("prior", "top_right"))
    )
  )
  expect_true(
    is_tibble(
      get_expected_category_statistic(
        fit, c("/d/", "/t/"), c("prior", "top_right"), c("mu", "Sigma")
      )
    )
  )
})

fit <- get_full_stanfit("NIW")

test_that("add ibbu draws - input check (3 cues)", {
  expect_true(is_tibble(get_draws(fit, groups = "prior")))
  expect_true(is_tibble(get_draws(fit, groups = "shifted_pos")))
  expect_true(is_tibble(get_draws(fit, groups = c("prior", "shifted_pos"))))
  expect_true(is_tibble(get_draws(fit, groups = "prior", ndraws = 1)))
  expect_error(get_draws(fit, groups = "priors"))
  expect_error(get_draws(fit, groups = "prior", ndraws = c(1, 2)))
})

test_that("add draws - output check (3 cues)", {
  expect_equal(
    length(
      unique(get_draws(fit, groups = "prior", ndraws = 10, seed = 1)$.draw)
    ),
    10
  )
  expect_equal(
    nrow(get_draws(fit, groups = "prior", summarize = TRUE)),
    length(get_category_labels(fit))
  )
  expect_equal(
    names(get_draws(fit, groups = "prior", summarize = TRUE)),
    c(".chain", ".iteration", ".draw", "group", "category", "kappa", "nu", "m",
      "S", "lapse_rate", "Sigma_exp", "Sigma_marg")
  )
  expect_equal(
    names(get_draws(fit, groups = "prior", summarize = TRUE, nest = TRUE)),
    c(".chain", ".iteration", ".draw", "group", "category", "kappa", "nu", "m",
      "S", "lapse_rate", "Sigma_exp", "Sigma_marg")
  )
  expect_equal(
    names(get_draws(fit, groups = "prior", summarize = TRUE, nest = FALSE)),
    c(".chain", ".iteration", ".draw", "group", "category", "kappa", "nu",
      "cue", "cue2", "m", "S", "lapse_rate", "Sigma_exp", "Sigma_marg")
  )
})

test_that("get exposure category statistic (NIW)", {
  expect_true(is.matrix(get_exposure_category_css(fit, "/d/", "no_exposure")))
  expect_error(
    is.matrix(get_exposure_category_css(fit, "wrong", "no_exposure"))
  )
  expect_error(is.matrix(get_exposure_category_css(fit, "/d/", "wrong")))
  expect_true(
    is_tibble(
      get_exposure_category_css(
        fit, c("/d/", "/t/"), c("no_exposure", "shifted_pos")
      )
    )
  )
  expect_true(is.matrix(get_exposure_category_uss(fit, "/d/", "no_exposure")))
  expect_error(
    is.matrix(get_exposure_category_uss(fit, "wrong", "no_exposure"))
  )
  expect_error(is.matrix(get_exposure_category_uss(fit, "/d/", "wrong")))
  expect_true(
    is_tibble(
      get_exposure_category_uss(
        fit, c("/d/", "/t/"), c("no_exposure", "shifted_pos")
      )
    )
  )
  expect_true(is.matrix(get_exposure_category_cov(fit, "/d/", "no_exposure")))
  expect_error(
    is.matrix(get_exposure_category_cov(fit, "wrong", "no_exposure"))
  )
  expect_error(is.matrix(get_exposure_category_cov(fit, "/d/", "wrong")))
  expect_true(
    is_tibble(
      get_exposure_category_cov(
        fit, c("/d/", "/t/"), c("no_exposure", "shifted_pos")
      )
    )
  )
  expect_true(
    is_tibble(
      get_exposure_category_statistic(
        fit,
        c("/d/", "/t/"),
        c("no_exposure", "shifted_pos"),
        c("n", "mean", "cov")
      )
    )
  )
})

test_that("get_parameter_names and get_params", {
  pars <- get_parameter_names(fit)
  expect_true(is.character(pars))
  expect_true(length(pars) > 0)
  expect_true("kappa_0" %in% pars || "m_0" %in% pars || "lp__" %in% pars)

  pars_orig <- get_parameter_names(fit, original_pars = TRUE)
  expect_true(is.character(pars_orig))
  expect_true(length(pars_orig) > 0)

  expect_warning(
    pars_dep <- get_params(fit),
    "deprecated"
  )
  expect_equal(pars_dep, pars)
})

test_that("get_draws warns when wide argument is supplied", {
  expect_warning(
    get_draws(fit, groups = "prior", ndraws = 2, seed = 1, wide = TRUE),
    class = "lifecycle_warning_deprecated"
  )
})
