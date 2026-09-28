test_that("plot_model_updates works for a list of cognitive models", {
  m1 <- example_mvg_ideal_observer()
  m2 <- example_mvg_ideal_observer()

  # 1D or 2D category update plot
  p_cat <- plot_model_updates(list("Initial" = m1, "Updated" = m2), what = "categories")
  expect_s3_class(p_cat, "ggplot")
  expect_true(".update_step" %in% names(p_cat$data))
  expect_equal(levels(p_cat$data$.update_step), c("Initial", "Updated"))

  # Categorization function update plot
  p_fn <- plot_model_updates(list("Initial" = m1, "Updated" = m2), what = "categorization_function")
  expect_s3_class(p_fn, "ggplot")
  expect_true(".update_step" %in% names(p_fn$data))
})

test_that("plot_model_updates works for MVBU_Stanfit objects", {
  fit <- get_example_stanfit(
    1,
    stanmodel = "NIW_ideal_adaptor"
  )

  p_stan_cat <- plot_model_updates(fit, what = "categories")
  expect_s3_class(p_stan_cat, "ggplot")
  expect_true(".update_step" %in% names(p_stan_cat$data))

  p_stan_fn <- plot_model_updates(fit, what = "categorization_function")
  expect_s3_class(p_stan_fn, "ggplot")
  expect_true(".update_step" %in% names(p_stan_fn$data))
})
