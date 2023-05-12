context("get information from NIW")

test_that("Get posterior predictive - input check x", {
  expect_error(get_NIW_posterior_predictive(1,
                                          matrix(c(0,0), nrow = 1),
                                          matrix(c(1, -.5, -.5, 2), nrow = 2),
                                          kappa = 1000,
                                          nu = 1000))
  expect_error(get_NIW_posterior_predictive(matrix(c(1,0,1,2,1,2), nrow = 2),
                                          matrix(c(0,0), nrow = 1),
                                          matrix(c(1, -.5, -.5, 2), nrow = 2),
                                          kappa = 1000,
                                          nu = 1000))
})


test_that("Get posterior predictive - input check S", {
  expect_error(get_NIW_posterior_predictive(matrix(c(1,0), nrow = 1),
                                        matrix(c(0,0), nrow = 1),
                                        matrix(c(1, 2, .5, 2), nrow = 2),
                                        kappa = 1000,
                                        nu = 1000),
               "sigma must be a symmetric matrix")
})


test_that("Get posterior predictive - output check", {
  expect_equal(1, length(get_NIW_posterior_predictive(x = matrix(c(1,0), nrow = 1),
                                                  m = matrix(c(0,0), nrow = 1),
                                                  S = matrix(c(1, -.5, -.5, 2), nrow = 2),
                                                  kappa = 1000,
                                                  nu = 1000)))
  expect_equal(2, length(get_NIW_posterior_predictive(x = matrix(c(1,0,2,-2), nrow = 2),
                                                  m = matrix(c(0,0), nrow = 1),
                                                  S = matrix(c(1, -.5, -.5, 2), nrow = 2),
                                                  kappa = 1000,
                                                  nu = 1000)))
  expect_equal(3, length(get_NIW_posterior_predictive(x = matrix(c(1,0,1,2,1,2), nrow = 3),
                                                  m = matrix(c(0,0), nrow = 1),
                                                  S = matrix(c(1, -.5, -.5, 2), nrow = 2),
                                                  kappa = 1000,
                                                  nu = 1000)))
})


