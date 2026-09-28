# Performance regression checks and cache/memory bloat safeguards
# Verifies copy-on-modify cache semantics, memory footprint limits,
# and computational reuse across cognitive models and Stanfit objects.

test_that("cognitive model closure cache maintains strict copy-on-modify isolation", {
  m1 <- example_model("MVG", n_cues = 2L)
  m2 <- example_model("MVG", n_cues = 2L)
  
  # 1. Fresh models have empty caches
  expect_equal(length(m1@cache), 0L)
  expect_equal(length(m2@cache), 0L)
  
  # 2. Extracting closure caches it on m1
  lf1 <- get_category_likelihood_function(m1)
  expect_true(is.function(lf1))
  
  # m2's cache must remain empty (no leakage across instances)
  expect_equal(length(m2@cache), 0L)
  
  # 3. Repeated retrieval returns the identical cached closure
  lf1_repeat <- get_category_likelihood_function(m1)
  expect_identical(lf1, lf1_repeat)
  
  pf1 <- get_category_posterior_function(m1)
  pf1_repeat <- get_category_posterior_function(m1)
  expect_identical(pf1, pf1_repeat)
})

test_that("cached cognitive model objects have bounded memory footprint (no bloat)", {
  m <- example_model("MVG", n_cues = 2L)
  
  # Populate likelihood and posterior closure cache
  lf <- get_category_likelihood_function(m)
  pf <- get_category_posterior_function(m)
  
  # Evaluate on test observations
  cues <- get_cue_labels(m)
  dat <- as.data.frame(matrix(rnorm(20L), nrow = 10L))
  colnames(dat) <- cues
  res_post <- pf(dat)
  res_crit <- categorize(m, dat, decision_rule = "criterion")
  
  # Cache slot should only contain known function closures, not raw evaluation datasets
  cache_names <- names(m@cache)
  expect_false(any(c("data", "observations", "grid", "evaluations") %in% cache_names))
  expect_true(all(vapply(m@cache, is.function, logical(1L))))
})

test_that("plotting functions do not leak grid surfaces into model cache", {
  m <- example_model("UVG", n_cues = 1L)
  
  # Render category plot and categorization function plot
  p_cat <- plot_categories(m, interactive = FALSE)
  p_dec <- plot_categorization_functions(m, interactive = FALSE)
  
  expect_s3_class(p_cat, "ggplot")
  expect_s3_class(p_dec, "ggplot")
  
  # Model cache must not retain large continuous grid evaluation matrices
  expect_false("grid" %in% names(m@cache))
  expect_false("surface" %in% names(m@cache))
})

test_that(
  "Stanfit caching reuses summaries and avoids redundant posterior extractions",
  {
    fit <- get_example_stanfit(
      1,
      stanmodel = "NIW_ideal_adaptor"
    )
  
  # 1. Add posterior latents caches draws and summary once
  fit_cached <- add_posterior_latents(fit)
  cache <- .get_cache(fit_cached)
  expect_false(is.null(cache$draws))
  expect_false(is.null(cache$summary))
  
  # 2. Getter calls reuse summary cleanly without altering structure
  summary_before <- cache$summary
  stat_exp <- get_expected_category_statistic(fit_cached, statistic = "mu")
  stat_marg <- get_marginal_category_statistic(fit_cached, statistic = "mu")
  
  expect_s3_class(stat_exp, "data.frame")
  expect_s3_class(stat_marg, "data.frame")
  expect_identical(.get_cache(fit_cached)$summary, summary_before)
})

test_that(
  "reconstruct_update_history produces bounded model lists efficiently",
  {
    fit <- get_example_stanfit(
      1,
      stanmodel = "NIW_ideal_adaptor"
    )
  
  # Fast mode (discard uncertainty)
  t0 <- proc.time()
  hist_fast <- reconstruct_update_history(
    fit,
    uncertainty_treatment = "discard",
    step_size = 20L
  )
  t_fast <- (proc.time() - t0)[["elapsed"]]
  
  expect_true(S7::S7_inherits(hist_fast, MVBU_ModelList))
  expect_true(length(hist_fast@models) >= 2L)
  expect_lt(t_fast, 5.0) # Must complete comfortably within 5 seconds
  
  # Each checkpoint model in fast mode is a valid CognitiveModel
  for (i in seq_along(hist_fast@models)) {
    expect_true(S7::S7_inherits(hist_fast[[i]], MVBU_CognitiveModel))
  }
})
