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

test_that("add_posterior_latents caches latents on Stanfit objects", {
  fit <- get_example_stanfit(
    1,
    stanmodel = "NIW_ideal_adaptor"
  )
  expect_null(.get_cache(fit)$draws)
  expect_null(.get_cache(fit)$summary)
  
  # Run post-processing helper
  fit_cached <- add_posterior_latents(fit)
  cache <- .get_cache(fit_cached)
  expect_false(is.null(cache$draws))
  expect_false(is.null(cache$summary))
  expect_s3_class(cache$draws$parameters, "tbl_df")
  expect_true(is.function(cache$draws$likelihood_function))
  expect_true(is.function(cache$draws$posterior_function))
})

test_that(
  "add_criterion and loo_compare work on Stanfit objects like brms",
  {
    fit1 <- get_example_stanfit(
      1,
      stanmodel = "NIW_ideal_adaptor"
    )
    fit1 <- add_criterion(
      fit1,
      criterion = c("loo", "waic", "bayes_R2"),
      model_name = "Model_NIW"
    )
  
  expect_true(all(c("loo", "waic", "bayes_R2") %in% names(fit1@criteria)))
  expect_s3_class(fit1@criteria$loo, "loo")
  expect_s3_class(fit1@criteria$waic, "waic")
  expect_s3_class(fit1@criteria$bayes_R2, "bayes_R2")
  expect_equal(fit1@metadata$model_name, "Model_NIW")

  # Test loo_compare on Stanfit objects
  fit2 <- get_example_stanfit(
    1,
    stanmodel = "NIW_ideal_adaptor"
  )
  comp <- loo_compare(fit1, fit2, criterion = "loo")
  expect_s3_class(comp, "compare.loo")
})

test_that("pp_check works on Stanfit objects", {
  fit <- get_example_stanfit(
    1,
    stanmodel = "NIW_ideal_adaptor"
  )
  skip_if_not_installed("bayesplot")
  
  p <- pp_check(fit, type = "dens_overlay", ndraws = 10)
  expect_s3_class(p, "ggplot")
})

