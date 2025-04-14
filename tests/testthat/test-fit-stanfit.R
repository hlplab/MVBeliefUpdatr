context("fit stanfit")

source("../functions-to-make-or-load-models.R")

test_that("Test fitting transformations (one cue)", {
  expect_no_error(fit <- get_example_stanfit(1, transform_type = "identity", control = list(adapt_delta = .95), silent = 0, verbose = T))
  expect_no_error(fit <- get_example_stanfit(1, transform_type = "center", control = list(adapt_delta = .95), silent = 0, verbose = T))
  expect_no_error(fit <- get_example_stanfit(1, transform_type = "standardize", control = list(adapt_delta = .95), silent = 0, verbose = T))
  expect_no_error(fit <- get_example_stanfit(1, transform_type = "PCA whiten", control = list(adapt_delta = .95), silent = 0, verbose = T))
  expect_no_error(fit <- get_example_stanfit(1, transform_type = "ZCA whiten", control = list(adapt_delta = .95), silent = 0, verbose = T))
  expect_error(fit <- get_example_stanfit(1, transform_type = "other", silent = 0, verbose = T))
})

test_that("Test fitting for more than one cue", {
  expect_no_error(fit <- get_example_stanfit(2, transform_type = "PCA whiten", control = list(adapt_delta = .975), silent = 0, verbose = T))
  expect_no_error(fit <- get_example_stanfit(3, transform_type = "PCA whiten", tau_scale = rep(2.5, 3), warmup = 1500, iter = 2500, control = list(adapt_delta = .995, max_treedepth = 15), silent = 0, verbose = T))
})

