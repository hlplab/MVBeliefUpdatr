test_that("plot_categories defaults to 2D slicing for 3-cue representations & models", {
  rep3d_A <- new_mvg_category_representation(
    category_labels = "A",
    cue_labels = c("F1", "F2", "F3"),
    mu = c(300, 1000, 2500),
    Sigma = diag(c(50^2, 100^2, 200^2))
  )
  rep3d_B <- new_mvg_category_representation(
    category_labels = "B",
    cue_labels = c("F1", "F2", "F3"),
    mu = c(500, 1500, 2700),
    Sigma = diag(c(60^2, 120^2, 220^2))
  )

  tpl3d <- new_category_representation_template(
    representations = list(A = rep3d_A, B = rep3d_B)
  )
  model3d <- new_mvg_ideal_observer(
    category_template = tpl3d,
    category_prior = c(A = 0.5, B = 0.5)
  )

  # 3-cue category plot defaults to 2D slicing
  p_cat <- plot_categories(model3d)
  expect_s3_class(p_cat, "ggplot")
  # Verify faceting exists
  expect_true(inherits(p_cat$facet, "FacetWrap"))

  # 3-cue categorization plot defaults to 2D slicing
  p_categorize <- plot_categorization_function(model3d)
  expect_s3_class(p_categorize, "ggplot")
  expect_true(inherits(p_categorize$facet, "FacetWrap"))

  # Custom slice cue and slice values
  p_custom <- plot_categories(
    model3d,
    slice_cue = "F1",
    slice_values = c(250, 400, 550)
  )
  expect_s3_class(p_custom, "ggplot")
})

test_that("parallel evaluation produces identical results to sequential evaluation", {
  mvg_rep1 <- new_mvg_category_representation(
    category_labels = "A",
    cue_labels = c("F1", "F2"),
    mu = c(0, 1),
    Sigma = diag(2)
  )
  mvg_rep2 <- new_mvg_category_representation(
    category_labels = "B",
    cue_labels = c("F1", "F2"),
    mu = c(2, 3),
    Sigma = diag(2)
  )
  tpl <- new_category_representation_template(
    representations = list(A = mvg_rep1, B = mvg_rep2)
  )
  model <- new_mvg_ideal_observer(
    category_template = tpl,
    category_prior = c(A = 0.5, B = 0.5)
  )

  grid <- expand.grid(
    F1 = seq(-2, 4, length.out = 30),
    F2 = seq(-1, 5, length.out = 30)
  )

  # Sequential vs parallel category densities
  dens_seq <- .mvbu_eval_category_densities(
    x = tpl,
    grid = grid,
    cues = c("F1", "F2"),
    parallel = FALSE
  )
  dens_par <- .mvbu_eval_category_densities(
    x = tpl,
    grid = grid,
    cues = c("F1", "F2"),
    parallel = TRUE,
    n_cores = 2L
  )
  expect_equal(dens_seq$density, dens_par$density, tolerance = 1e-9)

  # Sequential vs parallel categorization posteriors
  post_seq <- .mvbu_eval_categorization_posteriors(
    x = model,
    grid = grid,
    cues = c("F1", "F2"),
    parallel = FALSE
  )
  post_par <- .mvbu_eval_categorization_posteriors(
    x = model,
    grid = grid,
    cues = c("F1", "F2"),
    parallel = TRUE,
    n_cores = 2L
  )
  expect_equal(post_seq$posterior, post_par$posterior, tolerance = 1e-9)
})

test_that("plot_categorization and plot_correlations aliases work", {
  uvg_rep <- new_uvg_category_representation(
    category_labels = "A",
    cue_labels = "F1",
    mu = 0,
    sigma2 = 1
  )
  tpl <- new_category_representation_template(
    representations = list(A = uvg_rep)
  )
  model <- new_uvg_ideal_observer(
    category_template = tpl,
    category_prior = c(A = 1)
  )

  p_cat <- plot_categorization_function(model)
  expect_s3_class(p_cat, "ggplot")
})
