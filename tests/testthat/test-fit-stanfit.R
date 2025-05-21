source("../functions-to-make-or-load-models.R")

context("to_array")

x <- NULL
test_that("NULL inputs with simplify = TRUE", {
  expect_identical(
    to_array(x, simplify = TRUE), array(numeric()))
  expect_identical(
    to_array(x, inner_dims = 1, outer_dims = NULL, simplify = TRUE),
    array(numeric()))
  expect_identical(
    to_array(x, inner_dims = 1, outer_dims = 1, simplify = TRUE),
    array(numeric(), dim = c(0, 0)))
  expect_identical(
    to_array(x, inner_dims = 2, outer_dims = 1, simplify = TRUE),
    array(numeric(), dim = c(0, 0)))
  expect_identical(
    to_array(x, inner_dims = 2, outer_dims = 2, simplify = TRUE),
    array(numeric(), dim = c(0, 0)))
  expect_identical(
    to_array(x, inner_dims = c(1, 1), outer_dims = 1, simplify = TRUE),
    array(numeric(), dim = c(0, 0)))
  expect_identical(
    to_array(x, inner_dims = c(1, 1), outer_dims = c(1, 1), simplify = TRUE),
    array(numeric(), dim = c(0, 0, 0)))
})

test_that("NULL inputs with simplify = FALSE", {
  expect_identical(
    to_array(x, simplify = FALSE),
    array(numeric()))
  expect_identical(
    to_array(x, inner_dims = 1, outer_dims = NULL, simplify = FALSE),
    array(numeric()))
  expect_identical(
    to_array(x, inner_dims = 1, outer_dims = 1, simplify = FALSE),
    array(numeric(), dim = c(0, 0)))
  expect_identical(
    to_array(x, inner_dims = 2, outer_dims = 1, simplify = FALSE),
    array(numeric(), dim = c(0, 0)))
  expect_identical(
    to_array(x, inner_dims = 2, outer_dims = 2, simplify = FALSE),
    array(numeric(), dim = c(0, 0)))
  expect_identical(
    to_array(x, inner_dims = c(1, 1), outer_dims = 1, simplify = FALSE),
    array(numeric(), dim = c(0, 0, 0)))
  expect_identical(
    to_array(x, inner_dims = c(1, 1), outer_dims = c(1, 1), simplify = FALSE),
    array(numeric(), dim = c(0, 0, 0, 0)))
})

x <- 1
test_that("scalar inputs with simplify = TRUE", {
  expect_identical(
    to_array(x, simplify = TRUE), array(x))
  expect_identical(
    to_array(x, inner_dims = 1, outer_dims = NULL, simplify = TRUE),
    array(x))
  expect_error(
    to_array(x, inner_dims = 3, outer_dims = NULL, simplify = TRUE))
  expect_identical(
    to_array(x, inner_dims = 1, outer_dims = 3, simplify = TRUE),
    array(x, dim = c(3, 1)))
  expect_error(
    to_array(x, inner_dims = 2, outer_dims = 3, simplify = TRUE))
  expect_error(
    to_array(x, inner_dims = 2, outer_dims = 3, simplify = TRUE))
  # No simplification for atomic values
  expect_identical(
    to_array(x, inner_dims = c(1, 1), outer_dims = 3, simplify = TRUE),
    array(x, dim = c(3, 1, 1)))
  expect_identical(
    to_array(x, inner_dims = c(1, 1), outer_dims = c(3, 4), simplify = TRUE),
    array(x, dim = c(3, 4, 1, 1)))
})

test_that("scalar inputs with simplify = FALSE", {
  expect_identical(
    to_array(x, simplify = FALSE),
    array(x))
  expect_identical(
    to_array(x, inner_dims = 1, outer_dims = NULL, simplify = FALSE),
    array(x))
  expect_error(
    to_array(x, inner_dims = 3, outer_dims = NULL, simplify = FALSE))
  expect_identical(
    to_array(x, inner_dims = 1, outer_dims = 3, simplify = FALSE),
    array(x, dim = c(3, 1)))
  expect_error(
    to_array(x, inner_dims = 2, outer_dims = 3, simplify = FALSE))
  expect_error(
    to_array(x, inner_dims = 2, outer_dims = 3, simplify = FALSE))
  expect_identical(
    to_array(x, inner_dims = c(1, 1), outer_dims = 3, simplify = FALSE),
    array(x, dim = c(3, 1, 1)))
  expect_identical(
    to_array(x, inner_dims = c(1, 1), outer_dims = c(3, 4), simplify = FALSE),
    array(x, dim = c(3, 4, 1, 1)))
})

x <- 1:2
test_that("vector inputs with simplify = TRUE", {
  expect_identical(
    to_array(x, simplify = TRUE), array(x))
  expect_error(
    to_array(x, inner_dims = 1, outer_dims = NULL, simplify = TRUE))
  expect_identical(
    to_array(x, inner_dims = 2, outer_dims = NULL, simplify = TRUE),
    array(x))
  expect_equal(
    dim(to_array(x, inner_dims = 2, outer_dims = 3, simplify = TRUE)),
    c(3, 2))
  # No simplification for atomic values
  expect_error(
    to_array(x, inner_dims = c(1, 1), outer_dims = 3, simplify = TRUE))
  expect_error(
    to_array(x, inner_dims = c(1, 1), outer_dims = c(3, 4), simplify = TRUE))
})

test_that("vector inputs with simplify = FALSE", {
  expect_identical(
    to_array(x, simplify = FALSE), array(x))
  expect_error(
    to_array(x, inner_dims = 1, outer_dims = NULL, simplify = FALSE))
  expect_identical(
    to_array(x, inner_dims = 2, outer_dims = NULL, simplify = FALSE),
    array(x))
  expect_equal(
    dim(to_array(x, inner_dims = 2, outer_dims = 3, simplify = FALSE)),
    c(3, 2))
  expect_error(
    to_array(x, inner_dims = c(1, 1), outer_dims = 3, simplify = FALSE))
  expect_error(
    to_array(x, inner_dims = c(1, 1), outer_dims = c(3, 4), simplify = FALSE))
})

x <- matrix(1, nrow = 1)
test_that("matrix inputs with simplify = TRUE", {
  # No simplification for atomic values
  expect_equal(
    dim(to_array(x, simplify = TRUE)),
    c(1, 1))
  # Coercion
  expect_equal(
    dim(to_array(x, inner_dims = 1, outer_dims = NULL, simplify = TRUE)),
    c(1))
  expect_error(
    to_array(x, inner_dims = 2, outer_dims = NULL, simplify = TRUE))
  expect_error(
    to_array(x, inner_dims = 2, outer_dims = 3, simplify = TRUE))
  # No simplification for atomic values
  expect_equal(
    dim(to_array(x, inner_dims = c(1, 1), outer_dims = 3, simplify = TRUE)),
    c(3, 1, 1))
  expect_equal(
    dim(to_array(x, inner_dims = c(1, 1), outer_dims = c(3, 4), simplify = TRUE)),
    c(3, 4, 1, 1))
})

test_that("matrix inputs with simplify = FALSE", {
  # No simplification for atomic values
  expect_equal(
    dim(to_array(x, simplify = FALSE)),
    c(1, 1))
  # Coercion
  expect_equal(
    dim(to_array(x, inner_dims = 1, outer_dims = NULL, simplify = FALSE)),
    c(1))
  expect_error(
    to_array(x, inner_dims = 2, outer_dims = NULL, simplify = FALSE))
  expect_error(
    to_array(x, inner_dims = 2, outer_dims = 3, simplify = FALSE))
  expect_equal(
    dim(to_array(x, inner_dims = c(1, 1), outer_dims = 3, simplify = FALSE)),
    c(3, 1, 1))
  expect_equal(
    dim(to_array(x, inner_dims = c(1, 1), outer_dims = c(3, 4), simplify = FALSE)),
    c(3, 4, 1, 1))
})

# ADD LIST INPUT TEST


context("make_staninput_for_ideal_adaptor")

# Test whether input formatting works before testing fitting
test_that("NIX - test make_staninput_for_NIX_ideal_adaptor (one cue)", {
  expect_no_error(get_example_staninput(1, stanmodel = "NIX_ideal_adaptor", transform_type = "identity"))
  expect_error(get_example_staninput(1, stanmodel = "NIX_ideal_adaptor", transform_type = "center"))
  expect_error(get_example_staninput(1, stanmodel = "NIX_ideal_adaptor", transform_type = "standardize"))
  expect_error(get_example_staninput(1, stanmodel = "NIX_ideal_adaptor", transform_type = "PCA whiten"))
  expect_error(get_example_staninput(1, stanmodel = "NIX_ideal_adaptor", transform_type = "ZCA whiten"))
  expect_error(get_example_staninput(1, stanmodel = "NIX_ideal_adaptor", transform_type = "other"))
})

test_that("NIW - test make_staninput_for_NIW_ideal_adaptor (one cue)", {
  expect_no_error(get_example_staninput(1, stanmodel = "NIW_ideal_adaptor", transform_type = "identity"))
  expect_no_error(get_example_staninput(1, stanmodel = "NIW_ideal_adaptor", transform_type = "center"))
  expect_no_error(get_example_staninput(1, stanmodel = "NIW_ideal_adaptor", transform_type = "standardize"))
  expect_no_error(get_example_staninput(1, stanmodel = "NIW_ideal_adaptor", transform_type = "PCA whiten"))
  expect_no_error(get_example_staninput(1, stanmodel = "NIW_ideal_adaptor", transform_type = "ZCA whiten"))
  expect_error(get_example_staninput(1, stanmodel = "NIW_ideal_adaptor", transform_type = "other"))
})

test_that("MNIX - test make_staninput_for_MNIX_ideal_adaptor (one cue)", {
  expect_error(get_example_staninput(1, stanmodel = "MNIX_ideal_adaptor", transform_type = "identity"))
  expect_error(get_example_staninput(1, stanmodel = "MNIX_ideal_adaptor", transform_type = "center"))
  expect_error(get_example_staninput(1, stanmodel = "MNIX_ideal_adaptor", transform_type = "standardize"))
  expect_error(get_example_staninput(1, stanmodel = "MNIX_ideal_adaptor", transform_type = "PCA whiten"))
  expect_error(get_example_staninput(1, stanmodel = "MNIX_ideal_adaptor", transform_type = "ZCA whiten"))
  expect_error(get_example_staninput(1, stanmodel = "MNIX_ideal_adaptor", transform_type = "other"))
})

test_that("NIX - test make_staninput_for_NIX_ideal_adaptor (two cues)", {
  expect_error(get_example_staninput(2, stanmodel = "NIX_ideal_adaptor", transform_type = "identity"))
  expect_error(get_example_staninput(2, stanmodel = "NIX_ideal_adaptor", transform_type = "center"))
  expect_error(get_example_staninput(2, stanmodel = "NIX_ideal_adaptor", transform_type = "standardize"))
  expect_error(get_example_staninput(2, stanmodel = "NIX_ideal_adaptor", transform_type = "PCA whiten"))
  expect_error(get_example_staninput(2, stanmodel = "NIX_ideal_adaptor", transform_type = "ZCA whiten"))
  expect_error(get_example_staninput(2, stanmodel = "NIX_ideal_adaptor", transform_type = "other"))
})

test_that("NIW - test make_staninput_for_NIW_ideal_adaptor (two cues)", {
  expect_no_error(get_example_staninput(2, stanmodel = "NIW_ideal_adaptor", transform_type = "identity"))
  expect_no_error(get_example_staninput(2, stanmodel = "NIW_ideal_adaptor", transform_type = "center"))
  expect_no_error(get_example_staninput(2, stanmodel = "NIW_ideal_adaptor", transform_type = "standardize"))
  expect_no_error(get_example_staninput(2, stanmodel = "NIW_ideal_adaptor", transform_type = "PCA whiten"))
  expect_no_error(get_example_staninput(2, stanmodel = "NIW_ideal_adaptor", transform_type = "ZCA whiten"))
  expect_error(get_example_staninput(2, stanmodel = "NIW_ideal_adaptor", transform_type = "other"))
})

test_that("MNIX - test make_staninput_for_MNIX_ideal_adaptor (two cues)", {
  expect_no_error(get_example_staninput(2, stanmodel = "MNIX_ideal_adaptor", transform_type = "identity"))
  expect_no_error(get_example_staninput(2, stanmodel = "MNIX_ideal_adaptor", transform_type = "center"))
  expect_no_error(get_example_staninput(2, stanmodel = "MNIX_ideal_adaptor", transform_type = "standardize"))
  expect_no_error(get_example_staninput(2, stanmodel = "MNIX_ideal_adaptor", transform_type = "PCA whiten"))
  expect_no_error(get_example_staninput(2, stanmodel = "MNIX_ideal_adaptor", transform_type = "ZCA whiten"))
  expect_error(get_example_staninput(2, stanmodel = "MNIX_ideal_adaptor", transform_type = "other"))
})

# Test whether conditions with empty exposure information also work
test_that("NIX - test make_staninput_for_NIX_ideal_adaptor (two cues)", {
  expect_no_error(get_example_staninput(4, stanmodel = "NIX_ideal_adaptor", transform_type = "identity"))
  expect_error(get_example_staninput(4, stanmodel = "NIX_ideal_adaptor", transform_type = "ZCA whiten"))
})

test_that("NIW - test make_staninput_for_NIW_ideal_adaptor (two cues)", {
  expect_no_error(get_example_staninput(5, stanmodel = "NIW_ideal_adaptor", transform_type = "identity"))
  expect_no_error(get_example_staninput(5, stanmodel = "NIW_ideal_adaptor", transform_type = "ZCA whiten"))
})

test_that("MNIX - test make_staninput_for_MNIX_ideal_adaptor (two cues)", {
  expect_no_error(get_example_staninput(5, stanmodel = "MNIX_ideal_adaptor", transform_type = "identity"))
  expect_no_error(get_example_staninput(5, stanmodel = "MNIX_ideal_adaptor", transform_type = "ZCA whiten"))
})

# Test fitting
context("fit_ideal_adaptor")

test_that("NIX - test fitting transformations", {
  expect_no_error(fit <- get_example_stanfit(1, stanmodel = "NIX_ideal_adaptor", transform_type = "identity", control = list(adapt_delta = .95), silent = 0, verbose = F))
  # expect_error(fit <- get_example_stanfit(1, stanmodel = "NIX_ideal_adaptor", transform_type = "center", control = list(adapt_delta = .95), silent = 0, verbose = F))
  # expect_error(fit <- get_example_stanfit(1, stanmodel = "NIX_ideal_adaptor", transform_type = "standardize", control = list(adapt_delta = .95), silent = 0, verbose = F))
  # expect_error(fit <- get_example_stanfit(1, stanmodel = "NIX_ideal_adaptor", transform_type = "PCA whiten", control = list(adapt_delta = .95), silent = 0, verbose = F))
  # expect_error(fit <- get_example_stanfit(1, stanmodel = "NIX_ideal_adaptor", transform_type = "ZCA whiten", control = list(adapt_delta = .95), silent = 0, verbose = F))
  # expect_error(fit <- get_example_stanfit(1, stanmodel = "NIX_ideal_adaptor", transform_type = "other", silent = 0, verbose = T))
})

test_that("NIW - test fitting transformations (one cue)", {
  expect_no_error(fit <- get_example_stanfit(1, stanmodel = "NIW_ideal_adaptor", transform_type = "identity", control = list(adapt_delta = .95), silent = 0, verbose = F))
  expect_no_error(fit <- get_example_stanfit(1, stanmodel = "NIW_ideal_adaptor", transform_type = "center", control = list(adapt_delta = .95), silent = 0, verbose = F))
  expect_no_error(fit <- get_example_stanfit(1, stanmodel = "NIW_ideal_adaptor", transform_type = "standardize", control = list(adapt_delta = .95), silent = 0, verbose = F))
  expect_no_error(fit <- get_example_stanfit(1, stanmodel = "NIW_ideal_adaptor", transform_type = "PCA whiten", control = list(adapt_delta = .95), silent = 0, verbose = F))
  expect_no_error(fit <- get_example_stanfit(1, stanmodel = "NIW_ideal_adaptor", transform_type = "ZCA whiten", control = list(adapt_delta = .95), silent = 0, verbose = F))
  expect_error(fit <- get_example_stanfit(1, stanmodel = "NIW_ideal_adaptor", transform_type = "other", silent = 0, verbose = T))
})

test_that("NIX - test fitting for more than one cue", {
  expect_error(fit <- get_example_stanfit(2, stanmodel = "NIX_ideal_adaptor", transform_type = "PCA whiten", control = list(adapt_delta = .975), silent = 0, verbose = F))
  expect_error(fit <- get_example_stanfit(3, stanmodel = "NIX_ideal_adaptor", transform_type = "PCA whiten", tau_scale = rep(2.5, 3), warmup = 1500, iter = 2500,
                                             control = list(adapt_delta = .995, max_treedepth = 15), silent = 0, verbose = F))
})

test_that("NIW - test fitting for more than one cue", {
  expect_no_error(fit <- get_example_stanfit(2, stanmodel = "NIW_ideal_adaptor", transform_type = "PCA whiten", control = list(adapt_delta = .975), silent = 0, verbose = F))
  expect_no_error(fit <- get_example_stanfit(3, stanmodel = "NIW_ideal_adaptor", transform_type = "PCA whiten", tau_scale = rep(2.5, 3), warmup = 1500, iter = 2500,
                                             control = list(adapt_delta = .995, max_treedepth = 15), silent = 0, verbose = F))
})

test_that("MNIX - test fitting for more than one cue", {
  expect_no_error(fit <- get_example_stanfit(2, stanmodel = "MNIX_ideal_adaptor", transform_type = "PCA whiten", control = list(adapt_delta = .975), silent = 0, verbose = F))
  expect_no_error(fit <- get_example_stanfit(3, stanmodel = "MNIX_ideal_adaptor", transform_type = "PCA whiten", tau_scale = rep(2.5, 3), warmup = 1500, iter = 2500,
                                             control = list(adapt_delta = .995, max_treedepth = 15), silent = 0, verbose = F))
})

# TO DO: add testing for known mu_0 etc.
