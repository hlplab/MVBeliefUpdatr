context("get information from stanfit")

source("../functions-to-make-or-load-models.R")
fit1 <- get_example_stanfit(1)

test_that("Test is.NIW_ideal_adaptor_stanfit", {
  expect_false(is.NIW_ideal_adaptor_stanfit(NULL))
  expect_false(is.NIW_ideal_adaptor_stanfit(NA))
  expect_false(is.NIW_ideal_adaptor_stanfit(1))
  expect_false(is.NIW_ideal_adaptor_stanfit("1"))
  expect_false(is.NIW_ideal_adaptor_stanfit(TRUE))
  expect_false(is.NIW_ideal_adaptor_stanfit(list(1)))
  expect_false(is.NIW_ideal_adaptor_stanfit(example_exemplar_model(1)))
  expect_false(is.NIW_ideal_adaptor_stanfit(example_MVG_ideal_observer(1)))
  expect_false(is.NIW_ideal_adaptor_stanfit(example_NIW_ideal_adaptor(1)))
  expect_true(is.NIW_ideal_adaptor_stanfit(fit1))
})

test_that("add ibbu draws - input check", {
  expect_true(is_tibble(add_ibbu_stanfit_draws(fit1, groups = "prior")))
  expect_true(is_tibble(add_ibbu_stanfit_draws(fit1, groups = "plus2.2")))
  expect_true(is_tibble(add_ibbu_stanfit_draws(fit1, groups = c("prior", "plus2.2"))))
  expect_error(add_ibbu_stanfit_draws(fit1, groups = "priors"))
  expect_error(add_ibbu_stanfit_draws(fit1, groups = "prior", ndraws = 1))
  expect_error(add_ibbu_stanfit_draws(fit1, groups = "prior", ndraws = c(1, 2)))
})

test_that("add ibbu draws - output check", {
  expect_equal(nrow(add_ibbu_stanfit_draws(fit1, groups = "prior", ndraws = 10, seed = 1, wide = F) %>% distinct(.draw)), 10)
  expect_equal(nrow(add_ibbu_stanfit_draws(fit1, groups = "prior", wide = F, summarize = T)), 2)
  expect_equal(names(add_ibbu_stanfit_draws(fit1, groups = "prior", summarize = T)),
               c(".chain", ".iteration", ".draw", "group", "category", "kappa", "nu", "lapse_rate", "m", "S"))
  expect_equal(names(add_ibbu_stanfit_draws(fit1, groups = "prior", summarize = T, nest = T)),
               c(".chain", ".iteration", ".draw", "group", "category", "kappa", "nu", "lapse_rate", "m", "S"))
  expect_equal(names(add_ibbu_stanfit_draws(fit1, groups = "prior", summarize = T, nest = F)),
               c("cue", "cue2", ".chain", ".iteration", ".draw", "group", "category", "kappa", "nu", "lapse_rate", "m", "S"))
})

# test_that("Add ibbu draws - check wide = T", {
#   expect_equal(nrow(add_ibbu_stanfit_draws(fit1, groups = "prior", wide = T, summarize = T)), 1)
#   expect_equal(nrow(add_ibbu_stanfit_draws(fit1, groups = "prior", wide = T, summarize = F)), 10)
# })

test_that("get exposure category statistic", {
  # Get error when *exposure* statistics is requested for *prior*
  expect_error(get_exposure_category_statistic_from_stanfit(fit1, groups = "prior"))
  expect_error(get_exposure_mean_from_stanfit(fit1, groups = "prior"))
  expect_error(get_exposure_ss_from_stanfit(fit1, groups = "prior"))
  # When prior is not requested
  # get_exposure_mean
  expect_true(is.vector(get_exposure_mean_from_stanfit(fit1, "A", "baseline")))
  expect_error(is.vector(get_exposure_mean_from_stanfit(fit1, "wrong", "baseline")))
  expect_error(is.vector(get_exposure_mean_from_stanfit(fit1, "A", "wrong")))
  expect_true(is_tibble(get_exposure_mean_from_stanfit(fit1, c("A", "B"), c("baseline", "plus2.2"))))
  # get_exposure_css
  expect_true(is.matrix(get_exposure_css_from_stanfit(fit1, "A", "baseline")))
  expect_error(is.matrix(get_exposure_css_from_stanfit(fit1, "wrong", "baseline")))
  expect_error(is.matrix(get_exposure_css_from_stanfit(fit1, "A", "wrong")))
  expect_true(is_tibble(get_exposure_css_from_stanfit(fit1, c("A", "B"), c("baseline", "plus2.2"))))
  # get_exposure_uss
  expect_true(is.matrix(get_exposure_uss_from_stanfit(fit1, "A", "baseline")))
  expect_error(is.matrix(get_exposure_uss_from_stanfit(fit1, "wrong", "baseline")))
  expect_error(is.matrix(get_exposure_uss_from_stanfit(fit1, "A", "wrong")))
  expect_true(is_tibble(get_exposure_uss_from_stanfit(fit1, c("A", "B"), c("baseline", "plus2.2"))))
  # get_exposure_cov
  expect_true(is.matrix(get_exposure_cov_from_stanfit(fit1, "A", "baseline")))
  expect_error(is.matrix(get_exposure_cov_from_stanfit(fit1, "wrong", "baseline")))
  expect_error(is.matrix(get_exposure_cov_from_stanfit(fit1, "A", "wrong")))
  expect_true(is_tibble(get_exposure_cov_from_stanfit(fit1, c("A", "B"), c("baseline", "plus2.2"))))
  # get multiple exposure statistics
  expect_true(is_tibble(get_exposure_category_statistic_from_stanfit(fit1, c("A", "B"), c("baseline", "plus2.2"), c("n", "mean", "cov"))))
})

test_that("get expected category statistic", {
  expect_true(is.vector(get_expected_mu_from_stanfit(fit1, "A", "prior")))
  expect_error(is.vector(get_expected_mu_from_stanfit(fit1, "wrong", "prior")))
  expect_error(is.vector(get_expected_mu_from_stanfit(fit1, "A", "wrong")))
  expect_true(is.matrix(get_expected_sigma_from_stanfit(fit1, "A", "prior")))
  expect_error(is.matrix(get_expected_sigma_from_stanfit(fit1, "wrong", "prior")))
  expect_error(is.matrix(get_expected_sigma_from_stanfit(fit1, "A", "wrong")))
  expect_true(is_tibble(get_expected_sigma_from_stanfit(fit1, c("A", "B"), c("prior", "plus2.2"))))
  expect_true(is_tibble(get_expected_category_statistic_from_stanfit(fit1, c("A", "B"), c("prior", "plus2.2"), c("mu", "Sigma"))))
})





