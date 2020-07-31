context("Get information from stanfit")

fit = readRDS("../test models/IBBU_mv_fit_3 exposure groups_2 categories_2 cues_Drouin et al 2016.rds")

test_that("Test fit class", {
  expect_true(is.mv_ibbu_stanfit(fit))
})

test_that("Add ibbu draws - input check", {
  expect_failure(expect_error(add_ibbu_stanfit_draws(fit, which = "prior")))
  expect_failure(expect_error(add_ibbu_stanfit_draws(fit, which = "posterior")))
  expect_failure(expect_error(add_ibbu_stanfit_draws(fit, which = "both")))
  expect_error(add_ibbu_stanfit_draws(fit, which = "priors"))
  expect_error(add_ibbu_stanfit_draws(fit, which = "prior", draws = -1:1))
})

test_that("Add ibbu draws - check wide = T", {
  expect_equal(nrow(add_ibbu_stanfit_draws(fit, which = "prior", draws = 1:10, wide = T, summarize = T)), 1)
  expect_equal(nrow(add_ibbu_stanfit_draws(fit, which = "prior", draws = 1:10, wide = T, summarize = F)), 10)
})

test_that("Add ibbu draws - output check", {
  expect_equal(nrow(add_ibbu_stanfit_draws(fit, which = "prior", draws = 1:10, wide = F)), 20)
  expect_equal(nrow(add_ibbu_stanfit_draws(fit, which = "prior", draws = 1:10, wide = F, summarize = T)), 2)
  expect_equal(names(add_ibbu_stanfit_draws(fit, which = "prior", summarize = T)),
               c(".chain", ".iteration", ".draw", "group", "category", "kappa", "nu", "M", "S", "lapse_rate"))
  expect_equal(names(add_ibbu_stanfit_draws(fit, which = "prior", summarize = T, nest = T)),
               c(".chain", ".iteration", ".draw", "group", "category", "kappa", "nu", "M", "S", "lapse_rate"))
  expect_equal(names(add_ibbu_stanfit_draws(fit, which = "prior", summarize = T, nest = F)),
               c("cue", "cue2", ".chain", ".iteration", ".draw", "group", "category", "kappa", "nu", "M", "S", "lapse_rate"))
})

test_that("Get expected category statistic", {
  expect_success(get_expected_mu(fit, "sh", "prior"))
  expect_success(get_expected_sigma(fit, "sh", "prior"))
  expect_success(get_expected_sigma(fit, c("s","sh"), c("prior", "control")))
  expect_success(get_expected_sigma(fit, c("s","sh"), c("prior", "Control")))
  expect_success(get_expected_category_statistic(fit, c("s","sh"), c("prior", "Control"), c("mu", "Sigma")))
})





