context("plot information from stanfit")

source("../functions-to-make-or-load-models.R")

# 2D stanfit with all conditions having exposure
fit1 <- get_example_stanfit(1)
# 2D stanfit with a condition without exposure
fit2 <- get_example_stanfit(2)
# 1D stanfit with a condition without exposure
fit3 <- get_example_stanfit(3)

test_that("print is.NIW_ideal_adaptor_stanfit", {
  expect_no_error(print(fit1))
  expect_no_error(print(fit2))
  expect_no_error(print(fit3))
})

test_that("plot_ibbu_stanfit_parameters", {
  expect_no_error(plot_ibbu_stanfit_parameters(fit1))
  expect_no_error(plot_ibbu_stanfit_parameters(fit2))
  expect_no_error(plot_ibbu_stanfit_parameters(fit3))
})

test_that("plot_ibbu_stanfit_parameter_correlations", {
  expect_no_error(plot_ibbu_stanfit_parameter_correlations(fit1))
  expect_no_error(plot_ibbu_stanfit_parameter_correlations(fit2))
  expect_no_error(plot_ibbu_stanfit_parameter_correlations(fit3))
})

test_that("plot_expected_ibbu_stanfit_categories_density2D", {
  expect_no_error(plot_expected_ibbu_stanfit_categories_2D(fit1, type = "contour"))
  expect_no_error(plot_expected_ibbu_stanfit_categories_2D(fit1, type = "density"))
  expect_no_error(plot_expected_ibbu_stanfit_categories_2D(fit2, type = "contour"))
  expect_no_error(plot_expected_ibbu_stanfit_categories_2D(fit2, type = "density"))
})

test_that("plot_expected_ibbu_stanfit_categories_density2D", {
  expect_no_error(plot_expected_ibbu_stanfit_categories_contour2D(fit1))
  expect_no_error(plot_expected_ibbu_stanfit_categories_contour2D(fit2))
})

test_that("plot_expected_ibbu_stanfit_categories_density2D", {
  expect_no_error(plot_expected_ibbu_stanfit_categories_density2D(fit1))
  expect_no_error(plot_expected_ibbu_stanfit_categories_density2D(fit2))
})
