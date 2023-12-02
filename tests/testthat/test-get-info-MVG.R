#' @import curl
#' @import remotes

context("get information from MVG - INITIAL TEST ONLY")

library(curl)
if (has_internet()) remotes::install_github("joeystanley/joeysvowels")
library(joeysvowels)
data("idahoans")

my_model <- make_MVG_ideal_observer_from_data(idahoans, category = "vowel", cues = c("F1", "F2"))
x.1 <- idahoans %>% mutate(x = map(F1, ~ c(...))) %>% pull(x)
x.2 <- idahoans %>% mutate(x = map2(F1, F2, ~ c(...))) %>% pull(x)
x.3 <- idahoans %>% mutate(x = pmap(.l = list(F1, F2, F3), ~ c(...))) %>% pull(x)

test_that("Get MVG likelihood - input check x (single element, non-list)", {
  expect_error(get_MVG_likelihood(
    x = 1,
    mu = matrix(c(0,0), nrow = 1),
    Sigma = matrix(c(1, -.5, -.5, 2), nrow = 2)))
  expect_error(get_MVG_likelihood(
    x = matrix(c(1,0,1,2,1,2), nrow = 6),
    mu = matrix(c(0,0), nrow = 1),
    Sigma = matrix(c(1, -.5, -.5, 2), nrow = 2)))
  expect_error(get_MVG_likelihood(
    x = matrix(c(1,0,1,2,1,2), nrow = 2),
    mu = matrix(c(0,0), nrow = 1),
    Sigma = matrix(c(1, -.5, -.5, 2), nrow = 2)))
  expect_no_error(get_MVG_likelihood(
    x = matrix(c(1,0,1,2,1,2), nrow = 3),
    mu = matrix(c(0,0), nrow = 1),
    Sigma = matrix(c(1, -.5, -.5, 2), nrow = 2)))
})

test_that("Get MVG likelihood - input check x (single-element list)", {
  expect_error(get_MVG_likelihood(
    x = list(1),
    mu = matrix(c(0,0), nrow = 1),
    Sigma = matrix(c(1, -.5, -.5, 2), nrow = 2)))
  expect_error(get_MVG_likelihood(
    x = list(matrix(c(1,0,1,2,1,2), nrow = 6)),
    mu = matrix(c(0,0), nrow = 1),
    Sigma = matrix(c(1, -.5, -.5, 2), nrow = 2)))
  expect_error(get_MVG_likelihood(
    x = list(matrix(c(1,0,1,2,1,2), nrow = 2)),
    mu = matrix(c(0,0), nrow = 1),
    Sigma = matrix(c(1, -.5, -.5, 2), nrow = 2)))
  expect_no_error(get_MVG_likelihood(
    x = list(matrix(c(1,0,1,2,1,2), nrow = 3)),
    mu = matrix(c(0,0), nrow = 1),
    Sigma = matrix(c(1, -.5, -.5, 2), nrow = 2)))
})

test_that("Get MVG likelihood - input check x (multi-element list)", {
  expect_error(get_MVG_likelihood(
    x = map(rep(1, 3), ~ .x),
    mu = matrix(c(0,0), nrow = 1),
    Sigma = matrix(c(1, -.5, -.5, 2), nrow = 2)))
  expect_error(get_MVG_likelihood(
    x = list(matrix(c(1,0,1,2,1,2), nrow = 6), matrix(c(1,0,1,2,1,2), nrow = 6), matrix(c(1,0,1,2,1,2), nrow = 6)),
    mu = matrix(c(0,0), nrow = 1),
    Sigma = matrix(c(1, -.5, -.5, 2), nrow = 2)))
  expect_error(get_MVG_likelihood(
    x = list(matrix(c(1,0,1,2,1,2), nrow = 2), matrix(c(1,0,1,2,1,2), nrow = 2), matrix(c(1,0,1,2,1,2), nrow = 2)),
    mu = matrix(c(0,0), nrow = 1),
    Sigma = matrix(c(1, -.5, -.5, 2), nrow = 2)))
  expect_no_error(get_MVG_likelihood(
    x = list(matrix(c(1,0,1,2,1,2), nrow = 3), matrix(c(1,0,1,2,1,2), nrow = 3), matrix(c(1,0,1,2,1,2), nrow = 3)),
    mu = matrix(c(0,0), nrow = 1),
    Sigma = matrix(c(1, -.5, -.5, 2), nrow = 2)))
})

test_that("Get likelihood from MVG - input check x (single element, non-list)", {
  expect_error(get_likelihood_from_MVG(
    x = 1,
    model = my_model))
  expect_error(get_likelihood_from_MVG(
    x = matrix(c(1,0,1,2,1,2), nrow = 6),
    model = my_model))
  expect_error(get_likelihood_from_MVG(
    x = matrix(c(1,0,1,2,1,2), nrow = 2),
    model = my_model))
  expect_no_error(get_likelihood_from_MVG(
    x = matrix(c(1,0,1,2,1,2), nrow = 3),
    model = my_model))
})

test_that("Get likelihood from MVG - input check x (single-element list)", {
  expect_error(get_likelihood_from_MVG(
    x = list(1),
    model = my_model))
  expect_error(get_likelihood_from_MVG(
    x = list(matrix(c(1,0,1,2,1,2), nrow = 6)),
    model = my_model))
  expect_error(get_likelihood_from_MVG(
    x = list(matrix(c(1,0,1,2,1,2), nrow = 2)),
    model = my_model))
  expect_no_error(get_likelihood_from_MVG(
    x = list(matrix(c(1,0,1,2,1,2), nrow = 3)),
    model = my_model))
})

test_that("Get likelihood from MVG - input check x (multi-element list)", {
  expect_error(get_likelihood_from_MVG(
    x = map(rep(1, 3), ~ .x),
    model = my_model))
  expect_error(get_likelihood_from_MVG(
    x = list(matrix(c(1,0,1,2,1,2), nrow = 6), matrix(c(1,0,1,2,1,2), nrow = 6), matrix(c(1,0,1,2,1,2), nrow = 6)),
    model = my_model))
  expect_error(get_likelihood_from_MVG(
    x = list(matrix(c(1,0,1,2,1,2), nrow = 2), matrix(c(1,0,1,2,1,2), nrow = 2), matrix(c(1,0,1,2,1,2), nrow = 2)),
    model = my_model))
  expect_no_error(get_likelihood_from_MVG(
    x = list(matrix(c(1,0,1,2,1,2), nrow = 3), matrix(c(1,0,1,2,1,2), nrow = 3), matrix(c(1,0,1,2,1,2), nrow = 3)),
    model = my_model))
})

test_that("Get likelihood from MVG - input check x", {
  expect_error(get_likelihood_from_MVG(x = x.1, model = my_model))
  expect_error(get_likelihood_from_MVG(x = x.3, model = my_model))
  expect_no_error(get_likelihood_from_MVG(x = x.2, model = my_model))
})

test_that("Get categorization from MVG ideal observer - input check x", {
  expect_no_error(get_categorization_from_MVG_ideal_observer(x = x.2, model = my_model,
                                                             noise_treatment = "no_noise",
                                                             lapse_treatment = "no_lapses",
                                                             decision_rule = "sampling"))
})
