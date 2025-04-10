context("get information from stanfit")

source("../functions-to-make-or-load-models.R")

# 1D stanfit with all conditions having exposure
fit <- get_example_stanfit(1)

test_that("Test for single cue", {
  expect_no_error(summary(fit))
  expect_no_error(summary(fit, only_prior = TRUE))
  expect_no_error(summary(fit, include_transformed_pars = TRUE))
  expect_true(is.data.frame(summary(fit)))
})

# 3D stanfit with a condition without exposure
fit <- get_example_stanfit(3)

test_that("Test for multiple cues", {
  expect_no_error(summary(fit))
  expect_no_error(summary(fit, only_prior = TRUE))
  expect_no_error(summary(fit, include_transformed_pars = TRUE))
  expect_true(is.data.frame(summary(fit)))
})
