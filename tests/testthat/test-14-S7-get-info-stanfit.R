
skip_if_not(file.exists(example_stanfit_path(1, stanmodel = "NIW_ideal_adaptor")), "cached example stanfit not generated")
fit <- get_example_stanfit(1, stanmodel = "NIW_ideal_adaptor", file_refit = "never")


test_that("add draws - input check (1 cue)", {
  expect_true(is_tibble(get_draws(fit, groups = "prior")))
  expect_true(is_tibble(get_draws(fit, groups = "plus20")))
  expect_true(is_tibble(get_draws(fit, groups = c("prior", "plus20"))))
  expect_true(is_tibble(get_draws(fit, groups = "prior", ndraws = 1)))
  expect_error(get_draws(fit, groups = "priors"))
  expect_error(get_draws(fit, groups = "prior", ndraws = c(1, 2)))
})

test_that("add draws - output check (1 cue)", {
  expect_equal(length(unique(get_draws(fit, groups = "prior", ndraws = 10, seed = 1)$.draw)), 10)
  expect_equal(nrow(get_draws(fit, groups = "prior", summarize = T)), length(get_category_labels(fit)))
  expect_equal(names(get_draws(fit, groups = "prior", summarize = T)),
               c(".chain", ".iteration", ".draw", "group", "category", "kappa", "nu", "m", "S", "lapse_rate"))
  expect_equal(names(get_draws(fit, groups = "prior", summarize = T, nest = T)),
               c(".chain", ".iteration", ".draw", "group", "category", "kappa", "nu", "m", "S", "lapse_rate"))
  expect_equal(names(get_draws(fit, groups = "prior", summarize = T, nest = F)),
               c(".chain", ".iteration", ".draw", "group", "category", "kappa", "nu", "cue", "cue2", "m", "S", "lapse_rate"))
})

fit <- get_example_stanfit(2, stanmodel = "NIW_ideal_adaptor", transform_type = "standardize", file_refit = "never")
test_that("add draws - input check (2 cues)", {
  expect_true(is_tibble(get_draws(fit, groups = "prior")))
  expect_true(is_tibble(get_draws(fit, groups = "plus20.20")))
  expect_true(is_tibble(get_draws(fit, groups = c("prior", "plus20.20"))))
  expect_true(is_tibble(get_draws(fit, groups = "prior", ndraws = 1)))
  expect_error(get_draws(fit, groups = "priors"))
  expect_error(get_draws(fit, groups = "prior", ndraws = c(1, 2)))
})

test_that("add draws - output check (2 cues)", {
  expect_equal(length(unique(get_draws(fit, groups = "prior", ndraws = 10, seed = 1)$.draw)), 10)
  expect_equal(nrow(get_draws(fit, groups = "prior", summarize = T)), length(get_category_labels(fit)))
  expect_equal(names(get_draws(fit, groups = "prior", summarize = T)),
               c(".chain", ".iteration", ".draw", "group", "category", "kappa", "nu", "m", "S", "lapse_rate"))
  expect_equal(names(get_draws(fit, groups = "prior", summarize = T, nest = T)),
               c(".chain", ".iteration", ".draw", "group", "category", "kappa", "nu", "m", "S", "lapse_rate"))
  expect_equal(names(get_draws(fit, groups = "prior", summarize = T, nest = F)),
               c(".chain", ".iteration", ".draw", "group", "category", "kappa", "nu", "cue", "cue2", "m", "S", "lapse_rate"))
})

test_that("get exposure category statistic", {
  # Get error when *exposure* statistics is requested for *prior*
  expect_error(get_exposure_category_statistic(fit, groups = "prior"))
  expect_error(get_exposure_category_mean(fit, groups = "prior"))
  expect_error(get_exposure_category_cov(fit, groups = "prior"))
  # When prior is not requested
  # get_exposure_mean
  expect_true(is.vector(get_exposure_category_mean(fit, "/b/", "baseline")))
  expect_error(is.vector(get_exposure_category_mean(fit, "wrong", "baseline")))
  expect_error(is.vector(get_exposure_category_mean(fit, "/b/", "wrong")))
  expect_true(is_tibble(get_exposure_category_mean(fit, c("/b/", "/p/"), c("baseline", "plus20.20"))))
  # get_exposure_css
  expect_true(is.matrix(get_exposure_category_css(fit, "/b/", "baseline")))
  expect_error(is.matrix(get_exposure_category_css(fit, "wrong", "baseline")))
  expect_error(is.matrix(get_exposure_category_css(fit, "/b/", "wrong")))
  expect_true(is_tibble(get_exposure_category_css(fit, c("/b/", "/p/"), c("baseline", "plus20.20"))))
  # get_exposure_uss
  expect_true(is.matrix(get_exposure_category_uss(fit, "/b/", "baseline")))
  expect_error(is.matrix(get_exposure_category_uss(fit, "wrong", "baseline")))
  expect_error(is.matrix(get_exposure_category_uss(fit, "/b/", "wrong")))
  expect_true(is_tibble(get_exposure_category_uss(fit, c("/b/", "/p/"), c("baseline", "plus20.20"))))
  # get_exposure_cov
  expect_true(is.matrix(get_exposure_category_cov(fit, "/b/", "baseline")))
  expect_error(is.matrix(get_exposure_category_cov(fit, "wrong", "baseline")))
  expect_error(is.matrix(get_exposure_category_cov(fit, "/b/", "wrong")))
  expect_true(is_tibble(get_exposure_category_cov(fit, c("/b/", "/p/"), c("baseline", "plus20.20"))))
  # get multiple exposure statistics
  expect_true(is_tibble(get_exposure_category_statistic(fit, c("/b/", "/p/"), c("baseline", "plus20.20"), c("n", "mean", "cov"))))
})

test_that("get expected category statistic", {
  expect_true(is.vector(get_expected_mu(fit, "/b/", "prior")))
  expect_error(is.vector(get_expected_mu(fit, "wrong", "prior")))
  expect_error(is.vector(get_expected_mu(fit, "/b/", "wrong")))
  expect_true(is.matrix(get_expected_sigma(fit, "/b/", "prior")))
  expect_error(is.matrix(get_expected_sigma(fit, "wrong", "prior")))
  expect_error(is.matrix(get_expected_sigma(fit, "/b/", "wrong")))
  expect_true(is_tibble(get_expected_sigma(fit, c("/b/", "/p/"), c("prior", "plus20.20"))))
  expect_true(is_tibble(get_expected_category_statistic(fit, c("/b/", "/p/"), c("prior", "plus20.20"), c("mu", "Sigma"))))
})

test_that("get exposure category statistic", {
  # Get error when *exposure* statistics is requested for *prior*
  expect_error(get_exposure_category_statistic(fit, groups = "prior"))
  expect_error(get_exposure_category_mean(fit, groups = "prior"))
  expect_error(get_exposure_category_cov(fit, groups = "prior"))
  # When prior is not requested
  # get_exposure_mean
  expect_true(is.vector(get_exposure_category_mean(fit, "/b/", "baseline")))
  expect_error(is.vector(get_exposure_category_mean(fit, "wrong", "baseline")))
  expect_error(is.vector(get_exposure_category_mean(fit, "/b/", "wrong")))
  expect_true(is_tibble(get_exposure_category_mean(fit, c("/b/", "/p/"), c("baseline", "plus20.20"))))
  # get_exposure_css
  expect_true(is.matrix(get_exposure_category_css(fit, "/b/", "baseline")))
  expect_error(is.matrix(get_exposure_category_css(fit, "wrong", "baseline")))
  expect_error(is.matrix(get_exposure_category_css(fit, "/b/", "wrong")))
  expect_true(is_tibble(get_exposure_category_css(fit, c("/b/", "/p/"), c("baseline", "plus20.20"))))
  # get_exposure_uss
  expect_true(is.matrix(get_exposure_category_uss(fit, "/b/", "baseline")))
  expect_error(is.matrix(get_exposure_category_uss(fit, "wrong", "baseline")))
  expect_error(is.matrix(get_exposure_category_uss(fit, "/b/", "wrong")))
  expect_true(is_tibble(get_exposure_category_uss(fit, c("/b/", "/p/"), c("baseline", "plus20.20"))))
  # get_exposure_cov
  expect_true(is.matrix(get_exposure_category_cov(fit, "/b/", "baseline")))
  expect_error(is.matrix(get_exposure_category_cov(fit, "wrong", "baseline")))
  expect_error(is.matrix(get_exposure_category_cov(fit, "/b/", "wrong")))
  expect_true(is_tibble(get_exposure_category_cov(fit, c("/b/", "/p/"), c("baseline", "plus20.20"))))
  # get multiple exposure statistics
  expect_true(is_tibble(get_exposure_category_statistic(fit, c("/b/", "/p/"), c("baseline", "plus20.20"), c("n", "mean", "cov"))))
})

test_that("get expected category statistic", {
  expect_true(is.vector(get_expected_mu(fit, "/b/", "prior")))
  expect_error(is.vector(get_expected_mu(fit, "wrong", "prior")))
  expect_error(is.vector(get_expected_mu(fit, "/b/", "wrong")))
  expect_true(is.matrix(get_expected_sigma(fit, "/b/", "prior")))
  expect_error(is.matrix(get_expected_sigma(fit, "wrong", "prior")))
  expect_error(is.matrix(get_expected_sigma(fit, "/b/", "wrong")))
  expect_true(is_tibble(get_expected_sigma(fit, c("/b/", "/p/"), c("prior", "plus20.20"))))
  expect_true(is_tibble(get_expected_category_statistic(fit, c("/b/", "/p/"), c("prior", "plus20.20"), c("mu", "Sigma"))))
})

fit <- get_example_stanfit(3, stanmodel = "NIW_ideal_adaptor", transform_type = "standardize", file_refit = "never")
test_that("add ibbu draws - input check (3 cues)", {
  expect_true(is_tibble(get_draws(fit, groups = "prior")))
  expect_true(is_tibble(get_draws(fit, groups = "plus20.20.20")))
  expect_true(is_tibble(get_draws(fit, groups = c("prior", "plus20.20.20"))))
  expect_true(is_tibble(get_draws(fit, groups = "prior", ndraws = 1)))
  expect_error(get_draws(fit, groups = "priors"))
  expect_error(get_draws(fit, groups = "prior", ndraws = c(1, 2)))
})

test_that("add draws - output check (3 cues)", {
  expect_equal(length(unique(get_draws(fit, groups = "prior", ndraws = 10, seed = 1)$.draw)), 10)
  expect_equal(nrow(get_draws(fit, groups = "prior", summarize = T)), length(get_category_labels(fit)))
  expect_equal(names(get_draws(fit, groups = "prior", summarize = T)),
               c(".chain", ".iteration", ".draw", "group", "category", "kappa", "nu", "m", "S", "lapse_rate"))
  expect_equal(names(get_draws(fit, groups = "prior", summarize = T, nest = T)),
               c(".chain", ".iteration", ".draw", "group", "category", "kappa", "nu", "m", "S", "lapse_rate"))
  expect_equal(names(get_draws(fit, groups = "prior", summarize = T, nest = F)),
               c(".chain", ".iteration", ".draw", "group", "category", "kappa", "nu", "cue", "cue2", "m", "S", "lapse_rate"))
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

