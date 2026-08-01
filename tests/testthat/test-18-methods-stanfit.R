
# 1D stanfit with all conditions having exposure
fit <- readRDS("../benchmark-1cue.rds")

test_that("Test for single cue", {
  expect_no_error(suppressWarnings(summary(fit)))
  expect_no_error(suppressWarnings(summary(fit, only_prior = TRUE)))
  expect_no_error(suppressWarnings(summary(fit, include_transformed_pars = TRUE)))
  expect_true(is.list(suppressWarnings(summary(fit))) || is.matrix(suppressWarnings(summary(fit))) || is.table(suppressWarnings(summary(fit))) || is.data.frame(suppressWarnings(summary(fit))))
})

# 3D stanfit with a condition without exposure
fit <- readRDS("../benchmark-3cues.rds")

test_that("Test for multiple cues", {
  expect_no_error(suppressWarnings(summary(fit)))
  expect_no_error(suppressWarnings(summary(fit, only_prior = TRUE)))
  expect_no_error(suppressWarnings(summary(fit, include_transformed_pars = TRUE)))
  expect_true(is.list(suppressWarnings(summary(fit))) || is.matrix(suppressWarnings(summary(fit))) || is.table(suppressWarnings(summary(fit))) || is.data.frame(suppressWarnings(summary(fit))))
})
