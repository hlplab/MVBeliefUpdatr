context("plot information from stanfit")

source("../functions-to-make-or-load-models.R")

# 1D stanfit with all conditions having exposure
fit1 <- get_example_stanfit(1)
# 2D stanfit with a condition without exposure
fit2 <- get_example_stanfit(2)
# 3D stanfit with a condition without exposure
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

# TO DO: add more tests of selecting categories, groups, cues, and/or parameters to plot
test_that("plot_ibbu_stanfit_parameter_correlations", {
  expect_no_error(plot_ibbu_stanfit_parameter_correlations(fit1))
  expect_no_error(plot_ibbu_stanfit_parameter_correlations(fit2))
  expect_no_error(plot_ibbu_stanfit_parameter_correlations(fit3))
})

test_that("plot_expected_ibbu_stanfit_categories_density", {
  expect_warning(plot_expected_ibbu_stanfit_categories(fit1, type = "contour"))
  expect_no_error(plot_expected_ibbu_stanfit_categories(fit1, type = "density", resolution = 5, ndraws = 5))
  expect_no_error(plot_expected_ibbu_stanfit_categories(fit2, type = "contour"))
  expect_no_error(plot_expected_ibbu_stanfit_categories(fit2, type = "density", resolution = 5, ndraws = 5))
  expect_warning(plot_expected_ibbu_stanfit_categories(fit3, type = "contour"))
  expect_warning(plot_expected_ibbu_stanfit_categories(fit3, type = "density"))
})

test_that("plot_expected_ibbu_stanfit_categories_contour", {
  expect_error(plot_expected_ibbu_stanfit_categories_contour2D(fit1))
  expect_no_error(plot_expected_ibbu_stanfit_categories_contour2D(fit2))
  expect_error(plot_expected_ibbu_stanfit_categories_contour2D(fit3))
})

test_that("plot_expected_ibbu_stanfit_categories_density", {
  expect_no_error(plot_expected_ibbu_stanfit_categories_density1D(fit1, resolution = 5, ndraws = 5))
  expect_error(plot_expected_ibbu_stanfit_categories_density1D(fit2, resolution = 5, ndraws = 5))
  expect_error(plot_expected_ibbu_stanfit_categories_density1D(fit3))
  expect_error(plot_expected_ibbu_stanfit_categories_density2D(fit1, resolution = 5, ndraws = 5))
  expect_no_error(plot_expected_ibbu_stanfit_categories_density2D(fit2, resolution = 5, ndraws = 5))
  expect_error(plot_expected_ibbu_stanfit_categories_density2D(fit3))
})
