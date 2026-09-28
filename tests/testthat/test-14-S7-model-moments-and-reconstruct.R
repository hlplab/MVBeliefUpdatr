# =============================================================================
# Category Expected & Marginal Moments, Summary, and Reconstruction
# =============================================================================

test_that(
  "Expected & marginal moments match math definitions for NIW, NIX, MNIX",
  {
  # 1. NIW Category Representation
  cues <- c("cue1", "cue2")
  m_vec <- c(1.0, 2.0)
  S_mat <- matrix(c(2.0, 0.5, 0.5, 3.0), 2L, 2L)
  kappa_val <- 4.0
  nu_val <- 10.0
  D <- 2L

  niw_rep <- new_niw_category_representation(
    category_labels = "CatA",
    cue_labels = cues,
    m = m_vec,
    S = S_mat,
    kappa = kappa_val,
    nu = nu_val
  )

  expected_sigma_true <- S_mat / (nu_val - D - 1)
  marginal_sigma_true <- ((kappa_val + 1) / kappa_val) * expected_sigma_true

  expect_equal(as.numeric(get_expected_mu(niw_rep)), m_vec)
  expect_equal(as.numeric(get_marginal_mu(niw_rep)), m_vec)
  expect_equal(as.matrix(get_expected_sigma(niw_rep)), expected_sigma_true)
  expect_equal(as.matrix(get_marginal_sigma(niw_rep)), marginal_sigma_true)

  # Check marginal category statistic
  marg_stat <- get_marginal_category_statistic(
    niw_rep,
    statistic = c("mu", "Sigma")
  )
  expect_equal(as.numeric(marg_stat$mu), m_vec)
  expect_equal(as.matrix(marg_stat$Sigma), marginal_sigma_true)

  # 2. NIX Category Representation (1D)
  nix_rep <- new_nix_category_representation(
    category_labels = "Cat1",
    cue_labels = "x",
    m = 2.5,
    sigma2 = 1.8,
    kappa = 5.0,
    nu = 8.0
  )
  nix_s <- 1.8 * 8.0
  nix_exp_sig <- matrix(nix_s / (8.0 - 2.0), 1L, 1L)
  nix_marg_sig <- ((5.0 + 1.0) / 5.0) * nix_exp_sig

  expect_equal(as.numeric(get_expected_mu(nix_rep)), 2.5)
  expect_equal(as.numeric(get_marginal_mu(nix_rep)), 2.5)
  expect_equal(as.matrix(get_expected_sigma(nix_rep)), nix_exp_sig)
  expect_equal(as.matrix(get_marginal_sigma(nix_rep)), nix_marg_sig)

  # 3. MNIX Category Representation
  mnix_rep <- new_mnix_category_representation(
    category_labels = "CatM",
    cue_labels = cues,
    m = m_vec,
    sigma2 = c(1.5, 2.5),
    kappa = c(4.0, 4.0),
    nu = c(8.0, 8.0)
  )
  mnix_exp_diag <- (c(1.5, 2.5) * 8.0) / (8.0 - 2.0)
  mnix_marg_diag <- ((4.0 + 1.0) / 4.0) * mnix_exp_diag

  expect_equal(as.numeric(get_expected_mu(mnix_rep)), m_vec)
  expect_equal(as.numeric(get_marginal_mu(mnix_rep)), m_vec)
  expect_equal(diag(as.matrix(get_expected_sigma(mnix_rep))), mnix_exp_diag)
  expect_equal(diag(as.matrix(get_marginal_sigma(mnix_rep))), mnix_marg_diag)
})

test_that(
  "summary outputs expected and marginal moment tables for models & templates",
  {
  repA <- new_niw_category_representation(
    category_labels = "A",
    cue_labels = c("c1", "c2"),
    m = c(0, 0),
    S = diag(2, 2L),
    kappa = 2,
    nu = 6
  )
  repB <- new_niw_category_representation(
    category_labels = "B",
    cue_labels = c("c1", "c2"),
    m = c(2, 2),
    S = diag(2, 2L),
    kappa = 2,
    nu = 6
  )
  tpl <- new_category_representation_template(list(repA, repB))
  model <- new_niw_ideal_adaptor(category_template = tpl)

  # Summary of representation
  out_rep <- utils::capture.output(summary(repA))
  expect_true(any(grepl("Expected Moments", out_rep)))
  expect_true(any(grepl("Marginal Moments", out_rep)))

  # Summary of template
  out_tpl <- utils::capture.output(summary(tpl))
  expect_true(any(grepl("Expected.*Moments", out_tpl)))
  expect_true(any(grepl("Marginal.*Moments", out_tpl)))

  # Summary of model
  out_mod <- utils::capture.output(summary(model))
  expect_true(any(grepl("Expected.*Moments", out_mod)))
  expect_true(any(grepl("Marginal.*Moments", out_mod)))
})

test_that(
  "Stanfit/StanfitPosterior expected and marginal moments work",
  {
    sf_obj <- get_example_stanfit(
      1,
      stanmodel = "NIW_ideal_adaptor"
    )

  # get_draws computes Sigma_exp and Sigma_marg without redundant Sigma
  d_nest <- get_draws(sf_obj, nest = TRUE, summarize = FALSE)
  expect_true("Sigma_exp" %in% names(d_nest))
  expect_true("Sigma_marg" %in% names(d_nest))
  expect_false("Sigma" %in% names(d_nest))

  d_unnest <- get_draws(sf_obj, nest = FALSE, summarize = FALSE)
  expect_true("Sigma_exp" %in% names(d_unnest))
  expect_true("Sigma_marg" %in% names(d_unnest))
  expect_false("Sigma" %in% names(d_unnest))
  expect_true(all(c("cue", "cue2") %in% names(d_unnest)))

  # Expected vs marginal moments on Stanfit
  exp_stat <- get_expected_category_statistic(
    sf_obj,
    statistic = c("mu", "Sigma")
  )
  marg_stat <- get_marginal_category_statistic(
    sf_obj,
    statistic = c("mu", "Sigma")
  )
  expect_true(is.data.frame(exp_stat))
  expect_true(is.data.frame(marg_stat))
  expect_equal(nrow(exp_stat), nrow(marg_stat))

  # Check caching via add_posterior_latents
  sf_obj <- add_posterior_latents(sf_obj)
  expect_false(is.null(sf_obj@cache$draws))
  expect_false(is.null(sf_obj@cache$summary))
  expect_true(is.function(sf_obj@cache$draws$likelihood_function))
  expect_true(is.function(sf_obj@cache$draws$posterior_function))
  expect_true(is.function(sf_obj@cache$summary$likelihood_function))
  expect_true(is.function(sf_obj@cache$summary$posterior_function))

  # Test closure return types
  cues <- get_cue_labels(sf_obj)
  test_grid <- as.data.frame(matrix(c(0, 1), nrow = 2L, ncol = length(cues)))
  colnames(test_grid) <- cues
  lik_3d <- sf_obj@cache$draws$likelihood_function(test_grid, draws = 1:5)
  expect_equal(length(dim(lik_3d)), 3L)
  expect_equal(dim(lik_3d)[1L], 2L)
  expect_equal(dim(lik_3d)[3L], 5L)

  post_3d <- sf_obj@cache$draws$posterior_function(test_grid, draws = 1:5)
  expect_equal(length(dim(post_3d)), 3L)
  expect_equal(dim(post_3d)[1L], 2L)
  expect_equal(dim(post_3d)[3L], 5L)

  post_sum_2d <- sf_obj@cache$summary$posterior_function(test_grid)
  expect_equal(length(dim(post_sum_2d)), 2L)
  expect_equal(nrow(post_sum_2d), 2L)

  # summary on Stanfit contains expected and marginal moments, combined in print
  sum_sf <- summary(sf_obj)
  expect_true(S7::S7_inherits(sum_sf, Summary_MVBU_Stanfit))
  expect_true(is.data.frame(sum_sf@expected_moments))
  expect_true(is.data.frame(sum_sf@marginal_moments))
  expect_true(
    all(c("mean", "sd", "2.5%", "50%", "97.5%") %in%
      names(sum_sf@expected_moments))
  )
  expect_true(all(c("mu", "Sigma_exp") %in% sum_sf@expected_moments$Parameter))
  expect_true(
    all(c("mu", "Sigma_marg") %in% sum_sf@marginal_moments$Parameter)
  )

  # Custom probs in summary
  sum_custom <- summary(sf_obj, probs = c(0.10, 0.90))
  expect_true(
    all(c("mean", "sd", "10%", "90%") %in% names(sum_custom@expected_moments))
  )
  expect_true(
    all(c("mean", "sd", "10%", "90%") %in% names(sum_custom@marginal_moments))
  )

  out_sum <- utils::capture.output(print(sum_sf))
  expect_true(any(grepl("Category moments", out_sum, ignore.case = TRUE)))

  # MVBU_StanfitPosterior methods
  post_obj <- as_MVBU_stanfit_posterior(sf_obj)
  exp_mu_post <- get_expected_mu(post_obj)
  marg_mu_post <- get_marginal_mu(post_obj)
  expect_true(is.list(exp_mu_post) || is.numeric(exp_mu_post))
  expect_true(is.list(marg_mu_post) || is.numeric(marg_mu_post))

  out_post_sum <- utils::capture.output(summary(post_obj))
  expect_true(any(grepl("Category moments", out_post_sum, ignore.case = TRUE)))
})

test_that(
  "plot_categories and plot_categorization_functions uncertainty treatment",
  {
    sf_obj <- get_example_stanfit(
      1,
      stanmodel = "NIW_ideal_adaptor"
    )

  # plot_categories: uncertainty_treatment = "marginalize" vs "sample"
  p_cat_marg <- plot_categories(sf_obj, uncertainty_treatment = "marginalize")
  expect_true(inherits(p_cat_marg, "ggplot"))

  p_cat_sample <- plot_categories(
    sf_obj,
    uncertainty_treatment = "sample",
    ndraws = 5L
  )
  expect_true(inherits(p_cat_sample, "ggplot"))

  # plot_categorization_functions: "marginalize" vs "sample"
  p_func_marg <- plot_categorization_functions(
    sf_obj,
    uncertainty_treatment = "marginalize",
    ndraws = 5L,
    resolution = 20L
  )
  expect_true(inherits(p_func_marg, "ggplot"))

  p_func_sample <- plot_categorization_functions(
    sf_obj,
    uncertainty_treatment = "sample",
    ndraws = 5L,
    resolution = 20L
  )
  expect_true(inherits(p_func_sample, "ggplot"))
})

test_that(
  "reconstruct_update_history marginalize and discard work with summaries",
  {
    sf_obj <- get_example_stanfit(
      1,
      stanmodel = "NIW_ideal_adaptor"
    )

  # uncertainty_treatment = "marginalize" with 5 draws
  hist_full <- reconstruct_update_history(
    sf_obj,
    uncertainty_treatment = "marginalize",
    ndraws = 5L,
    step_size = 20L
  )
  expect_true(S7::S7_inherits(hist_full, MVBU_ModelList))
  expect_true(S7::S7_inherits(hist_full[[1]], MVBU_StanfitPosterior))

  # uncertainty_treatment = "discard"
  hist_fast <- reconstruct_update_history(
    sf_obj,
    uncertainty_treatment = "discard",
    step_size = 20L
  )
  expect_true(S7::S7_inherits(hist_fast, MVBU_ModelList))
  expect_true(S7::S7_inherits(hist_fast[[1]], MVBU_CognitiveModel))

  # summary on ModelList
  out_list_sum <- utils::capture.output(summary(hist_fast))
  expect_true(any(grepl("Model Moments Summary across Steps", out_list_sum)))
})
