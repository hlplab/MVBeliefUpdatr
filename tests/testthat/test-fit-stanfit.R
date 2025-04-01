context("fit stanfit")

source("../functions-to-make-or-load-models.R")

test_that("Test fitting stanfit", {
  # Testing default model storage
  expect_no_error(get_example_stanfit(1, center.observations = T, scale.observations = T, stanmodel = NULL))
  # With scaling
  expect_no_error(get_example_stanfit(1, center.observations = T, scale.observations = T))
  expect_no_error(get_example_stanfit(2, center.observations = T, scale.observations = T))
  expect_no_error(get_example_stanfit(3, center.observations = T, scale.observations = T))
  # Without scaling
  expect_no_error(get_example_stanfit(1, center.observations = T, scale.observations = F))
  expect_no_error(get_example_stanfit(2, center.observations = T, scale.observations = F))
  expect_no_error(get_example_stanfit(3, center.observations = T, scale.observations = F))
  # Without scaling, old stanmodel
  expect_no_error(get_example_stanfit(2, center.observations = T, scale.observations = F, stanmodel = "mvg_conj_sufficient_stats_lapse"))
  # check that forcing multivariate updating works even if there is only one cue
  expect_no_error(get_example_stanfit(1, center.observations = T, scale.observations = T, use_univariate_updating = F, iter = 100))
  # forcing univariate updating should throw an error if and only if there is more than 1 cue
  expect_no_error(get_example_stanfit(2, center.observations = T, scale.observations = T, use_univariate_updating = F, iter = 100))
})

