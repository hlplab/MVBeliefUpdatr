skip_if_not(
  all(vapply(1:3, function(ex) {
    file.exists(example_stanfit_path(ex, stanmodel = "NIW_ideal_adaptor"))
  }, logical(1))),
  "cached example stanfits not generated"
)

# 1D stanfit with all conditions having exposure
fit1 <- get_example_stanfit(
  1,
  stanmodel = "NIW_ideal_adaptor",
  file_refit = "never"
)
# 2D stanfit with a condition without exposure
fit2 <- get_example_stanfit(
  2,
  stanmodel = "NIW_ideal_adaptor",
  file_refit = "never"
)
# 3D stanfit with a condition without exposure
fit3 <- get_example_stanfit(
  3,
  stanmodel = "NIW_ideal_adaptor",
  file_refit = "never"
)

test_that("plot_parameters on Stanfit objects", {
  expect_no_error(plot_parameters(fit1))
  expect_no_error(plot_parameters(fit2))
  expect_no_error(plot_parameters(fit3))
})

test_that("plot_parameter_correlations on Stanfit objects", {
  expect_no_error(plot_parameter_correlations(fit1))
  expect_no_error(plot_parameter_correlations(fit2))
  expect_no_error(plot_parameter_correlations(fit3))
})

test_that("plot_parameters_pairwise on Stanfit objects", {
  expect_no_error(plot_parameters_pairwise(fit1))
  expect_no_error(plot_parameters_pairwise(fit2))
  p_pair_multi <- plot_parameters_pairwise(
    fit2,
    groups = get_group_labels(fit2, include_prior = TRUE)
  )
  expect_s3_class(p_pair_multi, "ggplot")
})

test_that("plot_categories on Stanfit objects", {
  # 1D stanfit
  p1 <- plot_categories(fit1, ndraws = 5, resolution = 10)
  expect_s3_class(p1, "ggplot")

  # 2D stanfit with contour and fill aesthetics
  p2_contour <- plot_categories(
    fit2,
    aes = "contour",
    ndraws = 5,
    show_exposure_data = TRUE,
    show_test_data = TRUE
  )
  expect_s3_class(p2_contour, "ggplot")

  p2_fill <- plot_categories(
    fit2,
    aes = "fill",
    resolution = 10,
    ndraws = 5
  )
  expect_s3_class(p2_fill, "ggplot")

  # 3D stanfit with cue projection
  cues_3d <- get_cue_labels(fit3)
  p3_proj1d <- plot_categories(fit3, cues = cues_3d[1], ndraws = 5)
  expect_s3_class(p3_proj1d, "ggplot")

  p3_proj2d <- plot_categories(
    fit3,
    cues = cues_3d[1:2],
    aes = "contour",
    ndraws = 5
  )
  expect_s3_class(p3_proj2d, "ggplot")
})

test_that("plot_categorization_function on Stanfit objects", {
  p_cat1 <- plot_categorization_function(fit1, ndraws = 5)
  expect_s3_class(p_cat1, "ggplot")

  p_cat2 <- plot_categorization_function(fit2, ndraws = 5)
  expect_s3_class(p_cat2, "ggplot")
})

test_that("plot_expected_categories deprecated wrapper works", {
  expect_warning(
    plot_expected_categories(fit1, type = "density", ndraws = 5),
    "was deprecated in MVBeliefUpdatr 0.1.0"
  )
  expect_warning(
    plot_expected_categories(fit2, type = "contour"),
    "was deprecated in MVBeliefUpdatr 0.1.0"
  )
})
