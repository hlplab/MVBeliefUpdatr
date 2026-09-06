test_that("closure caching works on MVBU_CognitiveModel for likelihood and posterior", {
  model <- example_mvg_ideal_observer()
  
  # Factory calls create and cache closures
  lf1 <- get_category_likelihood_function(model)
  lf2 <- get_category_likelihood_function(model)
  expect_identical(lf1, lf2)
  
  pf1 <- get_category_posterior_function(model)
  pf2 <- get_category_posterior_function(model)
  expect_identical(pf1, pf2)
  
  # Evaluated probabilities are numerically identical to direct calculation
  cues <- get_cue_labels(model)
  test_data <- matrix(rep(10, length(cues)), nrow = 1)
  colnames(test_data) <- cues
  
  post_cached <- pf1(test_data)
  post_direct <- posterior(model, test_data)
  expect_equal(post_cached, post_direct, tolerance = 1e-12)
})

test_that("add_parameter_draws caches pre-extracted draw matrices on Stanfit objects", {
  fit_file <- testthat::test_path("models", "example-stanfit-NIW_ideal_adaptor-1-standardize-42.rds")
  skip_if_not(file.exists(fit_file), "Fixture file not found")
  
  fit <- readRDS(fit_file)
  expect_null(.get_cache(fit)$draw_matrices)
  
  # Run optional post-processing helper
  fit_cached <- add_parameter_draws(fit)
  expect_false(is.null(.get_cache(fit_cached)$draw_matrices))
  expect_true(is.data.frame(.get_cache(fit_cached)$draw_matrices$draws_df))
  expect_equal(.get_cache(fit_cached)$draw_matrices$n_draws, 4000L)
})

test_that("add_criterion and loo_compare work on Stanfit objects following brms conventions", {
  fit_file <- testthat::test_path("models", "example-stanfit-NIW_ideal_adaptor-1-standardize-42.rds")
  skip_if_not(file.exists(fit_file), "Fixture file not found")
  
  fit1 <- readRDS(fit_file)
  fit1 <- add_criterion(fit1, criterion = c("loo", "waic", "bayes_R2"), model_name = "Model_NIW")
  
  expect_true(all(c("loo", "waic", "bayes_R2") %in% names(fit1@criteria)))
  expect_s3_class(fit1@criteria$loo, "loo")
  expect_s3_class(fit1@criteria$waic, "waic")
  expect_s3_class(fit1@criteria$bayes_R2, "bayes_R2")
  expect_equal(fit1@metadata$model_name, "Model_NIW")

  # Test loo_compare on Stanfit objects
  fit2 <- readRDS(fit_file)
  comp <- loo_compare(fit1, fit2, criterion = "loo")
  expect_s3_class(comp, "compare.loo")
})

test_that("pp_check works on Stanfit objects", {
  fit_file <- testthat::test_path("models", "example-stanfit-NIW_ideal_adaptor-1-standardize-42.rds")
  skip_if_not(file.exists(fit_file), "Fixture file not found")
  
  fit <- readRDS(fit_file)
  skip_if_not_installed("bayesplot")
  
  p <- pp_check(fit, type = "dens_overlay", ndraws = 10)
  expect_s3_class(p, "ggplot")
})

