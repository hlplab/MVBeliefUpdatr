context("plot information from stanfit")

source("../functions-to-make-or-load-models.R")

# 1D stanfit with all conditions having exposure
fit1 <- get_example_stanfit(1)
# 2D stanfit with a condition without exposure
fit2 <- get_example_stanfit(2)
# 3D stanfit with a condition without exposure
fit3 <- get_example_stanfit(3)

test_that("plot_parameters.NIW_ideal_adaptor_stanfit", {
  expect_no_error(plot_parameters(fit1))
  expect_no_error(plot_parameters(fit2))
  expect_no_error(plot_parameters(fit3))
  expect_no_error(plot_parameters.NIW_ideal_adaptor_stanfit(fit1))
  expect_no_error(plot_parameters.NIW_ideal_adaptor_stanfit(fit2))
  expect_no_error(plot_parameters.NIW_ideal_adaptor_stanfit(fit3))
})

# TO DO: add more tests of selecting categories, groups, cues, and/or parameters to plot
test_that("plot_parameter_correlations.NIW_ideal_adaptor_stanfit", {
  expect_no_error(plot_parameter_correlations(fit1))
  expect_no_error(plot_parameter_correlations(fit2))
  expect_no_error(plot_parameter_correlations(fit3))
  expect_no_error(plot_parameter_correlations.NIW_ideal_adaptor_stanfit(fit1))
  expect_no_error(plot_parameter_correlations.NIW_ideal_adaptor_stanfit(fit2))
  expect_no_error(plot_parameter_correlations.NIW_ideal_adaptor_stanfit(fit3))
})

test_that("plot_expected_categories.NIW_ideal_adaptor_stanfit", {
  expect_warning(plot_expected_categories(fit1, type = "contour"))
  expect_no_error(plot_expected_categories(fit1, type = "density", resolution = 5, ndraws = 5))
  expect_no_error(plot_expected_categories(fit2, type = "contour"))
  expect_no_error(plot_expected_categories(fit2, type = "density", resolution = 5, ndraws = 5))
  expect_warning(plot_expected_categories(fit3, type = "contour"))
  expect_warning(plot_expected_categories(fit3, type = "density"))
})

test_that("plot_expected_categories_contour.NIW_ideal_adaptor_stanfit", {
  expect_error(plot_expected_categories_contour2D.NIW_ideal_adaptor_stanfit(fit1))
  expect_no_error(plot_expected_categories_contour2D.NIW_ideal_adaptor_stanfit(fit2))
  expect_error(plot_expected_categories_contour2D.NIW_ideal_adaptor_stanfit(fit3))
})

test_that("plot_expected_categories_density.NIW_ideal_adaptor_stanfit", {
  expect_no_error(plot_expected_categories_density1D.NIW_ideal_adaptor_stanfit(fit1, resolution = 5, ndraws = 5))
  expect_error(plot_expected_categories_density1D.NIW_ideal_adaptor_stanfit(fit2, resolution = 5, ndraws = 5))
  expect_error(plot_expected_categories_density1D.NIW_ideal_adaptor_stanfit(fit3))
  expect_error(plot_expected_categories_density2D.NIW_ideal_adaptor_stanfit(fit1, resolution = 5, ndraws = 5))
  expect_no_error(plot_expected_categories_density2D.NIW_ideal_adaptor_stanfit(fit2, resolution = 5, ndraws = 5))
  expect_error(plot_expected_categories_density2D.NIW_ideal_adaptor_stanfit(fit3))
})
