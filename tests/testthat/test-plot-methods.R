test_that("plot_categories works on 1D representations and models", {
  uvg_rep <- new_uvg_category_representation(
    category_labels = "A",
    cue_labels = "F1",
    mu = 0,
    sigma2 = 1
  )
  p1 <- plot_categories(uvg_rep)
  expect_s3_class(p1, "ggplot")

  tpl <- new_category_representation_template(
    representations = list(A = uvg_rep)
  )
  p2 <- plot_categories(tpl)
  expect_s3_class(p2, "ggplot")

  model <- new_uvg_ideal_observer(
    category_template = tpl,
    category_prior = c(A = 1)
  )
  p3 <- plot_categories(model)
  expect_s3_class(p3, "ggplot")

  # Test NIX representation
  nix_rep <- new_nix_category_representation(
    category_labels = "A",
    cue_labels = "F1",
    m = 0,
    sigma2 = 1,
    kappa = 10,
    nu = 15
  )
  p4 <- plot_categories(nix_rep)
  expect_s3_class(p4, "ggplot")
})

test_that("plot_categories works on 2D representations and models", {
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

  p1 <- plot_categories(mvg_rep1, aes = "contour")
  expect_s3_class(p1, "ggplot")

  p1_fill <- plot_categories(mvg_rep1, aes = "fill")
  expect_s3_class(p1_fill, "ggplot")

  tpl <- new_category_representation_template(
    representations = list(A = mvg_rep1, B = mvg_rep2)
  )
  p2 <- plot_categories(tpl, aes = "contour")
  expect_s3_class(p2, "ggplot")

  p2_fill <- plot_categories(tpl, aes = "fill")
  expect_s3_class(p2_fill, "ggplot")

  model <- new_mvg_ideal_observer(
    category_template = tpl,
    category_prior = c(A = 0.5, B = 0.5)
  )
  p3 <- plot_categories(model)
  expect_s3_class(p3, "ggplot")
})

test_that("plot_categories handles marginalization for D >= 3 cues", {
  rep3d <- new_mvg_category_representation(
    category_labels = "A",
    cue_labels = c("F1", "F2", "F3"),
    mu = c(0, 1, 2),
    Sigma = diag(3)
  )
  tpl3d <- new_category_representation_template(
    representations = list(A = rep3d)
  )
  model3d <- new_mvg_ideal_observer(
    category_template = tpl3d,
    category_prior = c(A = 1)
  )

  # Marginalize to 1D
  p_1d <- plot_categories(model3d, cues = "F2")
  expect_s3_class(p_1d, "ggplot")

  # Marginalize to 2D
  p_2d <- plot_categories(model3d, cues = c("F1", "F3"))
  expect_s3_class(p_2d, "ggplot")

  # Error on invalid cue
  expect_error(plot_categories(model3d, cues = "F99"))
})

test_that("plot_categorization_function works on 1D and 2D models", {
  mvg_rep1 <- new_mvg_category_representation(
    category_labels = "A",
    cue_labels = c("F1", "F2"),
    mu = c(0, 0),
    Sigma = diag(2)
  )
  mvg_rep2 <- new_mvg_category_representation(
    category_labels = "B",
    cue_labels = c("F1", "F2"),
    mu = c(2, 2),
    Sigma = diag(2)
  )
  tpl <- new_category_representation_template(
    representations = list(A = mvg_rep1, B = mvg_rep2)
  )
  model <- new_mvg_ideal_observer(
    category_template = tpl,
    category_prior = c(A = 0.5, B = 0.5)
  )

  # 1D categorization function
  p_cat1d <- plot_categorization_function(model, cues = "F1")
  expect_s3_class(p_cat1d, "ggplot")

  # 2D categorization function with contour and fill
  p_cat2d_contour <- plot_categorization_function(
    model,
    cues = c("F1", "F2"),
    aes = "contour"
  )
  expect_s3_class(p_cat2d_contour, "ggplot")

  p_cat2d_fill <- plot_categorization_function(
    model,
    cues = c("F1", "F2"),
    aes = "fill"
  )
  expect_s3_class(p_cat2d_fill, "ggplot")
})

test_that("plot_parameters works on cognitive models", {
  mvg_rep <- new_mvg_category_representation(
    category_labels = "A",
    cue_labels = c("F1", "F2"),
    mu = c(0, 1),
    Sigma = diag(2)
  )
  tpl <- new_category_representation_template(
    representations = list(A = mvg_rep)
  )
  model <- new_mvg_ideal_observer(
    category_template = tpl,
    category_prior = c(A = 1)
  )

  p <- plot_parameters(model)
  expect_s3_class(p, "ggplot")
})

test_that("base plot() dispatches to plot_categories", {
  uvg_rep <- new_uvg_category_representation(
    category_labels = "A",
    cue_labels = "F1",
    mu = 0,
    sigma2 = 1
  )
  expect_s3_class(plot(uvg_rep), "ggplot")

  tpl <- new_category_representation_template(
    representations = list(A = uvg_rep)
  )
  expect_s3_class(plot(tpl), "ggplot")

  model <- new_uvg_ideal_observer(
    category_template = tpl,
    category_prior = c(A = 1)
  )
  expect_s3_class(plot(model), "ggplot")
})

test_that("plot methods work on MVBU_Stanfit", {
  fit_file <- testthat::test_path(
    "models",
    "example-stanfit-NIW_ideal_adaptor-1-standardize-42.rds"
  )
  skip_if_not(file.exists(fit_file), "Fixture file not found")
  fit <- readRDS(fit_file)

  p_cat <- plot_categories(fit)
  expect_s3_class(p_cat, "ggplot")

  p_cat_data <- plot_categories(
    fit,
    show_exposure_data = TRUE,
    show_test_data = TRUE
  )
  expect_s3_class(p_cat_data, "ggplot")

  p_par <- plot_parameters(fit)
  expect_s3_class(p_par, "ggplot")

  p_par_list <- plot_parameters(fit, combine_into_single_plot = FALSE)
  expect_type(p_par_list, "list")

  p_diag <- plot_diagnostics(fit)
  expect_s3_class(p_diag, "ggplot")

  p_base <- plot(fit)
  expect_s3_class(p_base, "ggplot")
})

test_that("category subsetting and single target default work", {
  repA <- new_uvg_category_representation("A", "F1", mu = 300, sigma2 = 50^2)
  repB <- new_uvg_category_representation("B", "F1", mu = 600, sigma2 = 60^2)
  tpl <- new_category_representation_template(list(A = repA, B = repB))
  mod <- new_uvg_ideal_observer(tpl, category_prior = c(A = 0.5, B = 0.5))

  # plot_categories with category subset
  p_sub <- plot_categories(tpl, categories = "A")
  expect_s3_class(p_sub, "ggplot")

  # plot_categorization_function defaults to first category
  p_catfun_def <- plot_categorization_function(mod, cues = "F1")
  expect_s3_class(p_catfun_def, "ggplot")
  expect_identical(p_catfun_def$labels$y, "P(A | F1)")

  # explicit category selection
  p_catfun_b <- plot_categorization_function(
    mod,
    cues = "F1",
    categories = "B"
  )
  expect_s3_class(p_catfun_b, "ggplot")
  expect_identical(p_catfun_b$labels$y, "P(B | F1)")
})

test_that("levels list and exemplar sampling overlay work", {
  set.seed(42)
  ex_pts <- matrix(rnorm(300), ncol = 2)
  colnames(ex_pts) <- c("cue1", "cue2")
  ex_rep <- new_exemplar_category_representation(
    category_labels = "E",
    cue_labels = c("cue1", "cue2"),
    exemplars = ex_pts
  )
  tpl <- new_category_representation_template(list(E = ex_rep))
  mod <- new_exemplar_model(tpl, category_prior = c(E = 1))

  # list levels
  p_lvl <- plot_categories(
    tpl,
    aes = c("fill", "contour"),
    levels = list(contour = c(0.5, 0.95), fill = c(0.5, 0.8)),
    n_exemplars = 50
  )
  expect_s3_class(p_lvl, "ggplot")

  # multi-panel parameters with combine_into_single_plot = FALSE
  p_list <- plot_parameters(mod, combine_into_single_plot = FALSE)
  expect_type(p_list, "list")

  # multi-panel parameters combined
  p_comb <- plot_parameters(mod, combine_into_single_plot = TRUE)
  expect_s3_class(p_comb, "ggplot")
})

test_that(
  "limits argument works for plot_categories & categorization_function",
  {
  repA <- new_uvg_category_representation("A", "F1", mu = 300, sigma2 = 50^2)
  tpl1d <- new_category_representation_template(list(A = repA))
  mod1d <- new_uvg_ideal_observer(tpl1d, category_prior = c(A = 1))

  # 1D limits
  p_lim1d <- plot_categories(mod1d, limits = c(200, 400))
  expect_s3_class(p_lim1d, "ggplot")
  p_cat_lim1d <- plot_categorization_function(mod1d, limits = c(200, 400))
  expect_s3_class(p_cat_lim1d, "ggplot")

  # 2D limits
  mvg_rep1 <- new_mvg_category_representation(
    category_labels = "A",
    cue_labels = c("F1", "F2"),
    mu = c(0, 0),
    Sigma = diag(2)
  )
  tpl2d <- new_category_representation_template(list(A = mvg_rep1))
  mod2d <- new_mvg_ideal_observer(tpl2d, category_prior = c(A = 1))

  p_lim2d <- plot_categories(
    mod2d,
    limits = list(F1 = c(-2, 2), F2 = c(-3, 3))
  )
  expect_s3_class(p_lim2d, "ggplot")

  p_cat_lim2d <- plot_categorization_function(
    mod2d,
    limits = list(F1 = c(-2, 2), F2 = c(-3, 3))
  )
  expect_s3_class(p_cat_lim2d, "ggplot")
})

test_that("pairwise parameter distributions work on MVBU_Stanfit", {
  fit_file <- testthat::test_path(
    "models",
    "example-stanfit-NIW_ideal_adaptor-1-standardize-42.rds"
  )
  skip_if_not(file.exists(fit_file), "Fixture file not found")
  fit <- readRDS(fit_file)

  p_pair <- plot_parameters_pairwise(fit)
  expect_s3_class(p_pair, "ggplot")
  expect_identical(p_pair$labels$title, "Pairwise parameter distributions")

  p_corr <- plot_parameter_correlations(fit)
  expect_s3_class(p_corr, "ggplot")
})

test_that("deprecated plotting functions throw lifecycle warnings", {
  repA <- new_uvg_category_representation("A", "F1", mu = 300, sigma2 = 50^2)
  tpl <- new_category_representation_template(list(A = repA))
  mod <- new_uvg_ideal_observer(tpl, category_prior = c(A = 1))

  expect_warning(
    p_dep1 <- plot_expected_categories(mod),
    class = "lifecycle_warning_deprecated"
  )
  expect_s3_class(p_dep1, "ggplot")

  expect_warning(
    p_dep2 <- plot_expected_categorization_function_1D(mod),
    class = "lifecycle_warning_deprecated"
  )
  expect_s3_class(p_dep2, "ggplot")
})

