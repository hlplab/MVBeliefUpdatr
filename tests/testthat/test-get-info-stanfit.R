context("get information from stanfit")

source("../functions-to-make-or-load-models.R")
fit <- get_example_stanfit()

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
  expect_true(is.NIW_ideal_adaptor_stanfit(fit))
})

test_that("add ibbu draws - input check", {
  expect_true(is_tibble(add_ibbu_stanfit_draws(fit, groups = "prior")))
  expect_true(is_tibble(add_ibbu_stanfit_draws(fit, groups = "plus2.2")))
  expect_true(is_tibble(add_ibbu_stanfit_draws(fit, groups = c("prior", "plus2.2"))))
  expect_error(add_ibbu_stanfit_draws(fit, groups = "priors"))
  expect_error(add_ibbu_stanfit_draws(fit, groups = "prior", ndraws = 1))
  expect_error(add_ibbu_stanfit_draws(fit, groups = "prior", ndraws = c(1, 2)))
})

test_that("add ibbu draws - output check", {
  expect_equal(nrow(add_ibbu_stanfit_draws(fit, groups = "prior", ndraws = 10, seed = 1, wide = F) %>% distinct(.draw)), 10)
  expect_equal(nrow(add_ibbu_stanfit_draws(fit, groups = "prior", wide = F, summarize = T)), 2)
  expect_equal(names(add_ibbu_stanfit_draws(fit, groups = "prior", summarize = T)),
               c(".chain", ".iteration", ".draw", "group", "category", "kappa", "nu", "lapse_rate", "m", "S"))
  expect_equal(names(add_ibbu_stanfit_draws(fit, groups = "prior", summarize = T, nest = T)),
               c(".chain", ".iteration", ".draw", "group", "category", "kappa", "nu", "lapse_rate", "m", "S"))
  expect_equal(names(add_ibbu_stanfit_draws(fit, groups = "prior", summarize = T, nest = F)),
               c("cue", "cue2", ".chain", ".iteration", ".draw", "group", "category", "kappa", "nu", "lapse_rate", "m", "S"))
})

# test_that("Add ibbu draws - check wide = T", {
#   expect_equal(nrow(add_ibbu_stanfit_draws(fit, groups = "prior", wide = T, summarize = T)), 1)
#   expect_equal(nrow(add_ibbu_stanfit_draws(fit, groups = "prior", wide = T, summarize = F)), 10)
# })


test_that("get expected category statistic", {
  expect_true(is.vector(get_expected_mu_from_stanfit(fit, "A", "prior")))
  expect_error(is.vector(get_expected_mu_from_stanfit(fit, "wrong", "prior")))
  expect_error(is.vector(get_expected_mu_from_stanfit(fit, "A", "wrong")))
  expect_true(is.matrix(get_expected_sigma_from_stanfit(fit, "A", "prior")))
  expect_error(is.matrix(get_expected_sigma_from_stanfit(fit, "wrong", "prior")))
  expect_error(is.matrix(get_expected_sigma_from_stanfit(fit, "A", "wrong")))
  expect_true(is_tibble(get_expected_sigma_from_stanfit(fit, c("A","B"), c("prior", "plus2.2"))))
  expect_true(is_tibble(get_expected_category_statistic_from_stanfit(fit, c("A","B"), c("prior", "plus2.2"), c("mu", "Sigma"))))
})





