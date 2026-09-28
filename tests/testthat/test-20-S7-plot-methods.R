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

test_that("plot_categorization_functions works on 1D and 2D models", {
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
  p_cat1d <- plot_categorization_functions(model, cues = "F1")
  expect_s3_class(p_cat1d, "ggplot")

  # 2D categorization function with contour and fill
  p_cat2d_contour <- plot_categorization_functions(
    model,
    cues = c("F1", "F2"),
    aes = "contour"
  )
  expect_s3_class(p_cat2d_contour, "ggplot")

  p_cat2d_fill <- plot_categorization_functions(
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
  fit <- get_example_stanfit(
    1,
    stanmodel = "NIW_ideal_adaptor"
  )

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

  # plot_categorization_functions defaults to first category
  p_catfun_def <- plot_categorization_functions(mod, cues = "F1")
  expect_s3_class(p_catfun_def, "ggplot")
  expect_identical(p_catfun_def$labels$y, "P(A | F1)")

  # explicit category selection
  p_catfun_b <- plot_categorization_functions(
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
  p_cat_lim1d <- plot_categorization_functions(mod1d, limits = c(200, 400))
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

  p_cat_lim2d <- plot_categorization_functions(
    mod2d,
    limits = list(F1 = c(-2, 2), F2 = c(-3, 3))
  )
  expect_s3_class(p_cat_lim2d, "ggplot")
})

test_that("pairwise parameter distributions work on MVBU_Stanfit", {
  fit <- get_example_stanfit(
    1,
    stanmodel = "NIW_ideal_adaptor"
  )

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

test_that("plot_categories defaults to line-only contour and handles exemplar sampling subtitle", {
  repA <- new_uvg_category_representation("A", "F1", mu = 300, sigma2 = 50^2)
  tpl1d <- new_category_representation_template(list(A = repA))

  # Single rep and template both default to line-only (no ribbon layer)
  p_rep <- plot_categories(repA)
  expect_s3_class(p_rep, "ggplot")
  layer_types_rep <- vapply(p_rep$layers, function(l) class(l$geom)[1], character(1))
  expect_true("GeomLine" %in% layer_types_rep)
  expect_false("GeomRibbon" %in% layer_types_rep)

  p_tpl <- plot_categories(tpl1d)
  expect_s3_class(p_tpl, "ggplot")
  layer_types_tpl <- vapply(p_tpl$layers, function(l) class(l$geom)[1], character(1))
  expect_true("GeomLine" %in% layer_types_tpl)
  expect_false("GeomRibbon" %in% layer_types_tpl)

  # Template subtitle formatting states representation type without category count
  expect_match(deparse1(p_tpl$labels$subtitle), "1.*D.*univariate.*Gaussian.*category")
  expect_no_match(deparse1(p_tpl$labels$subtitle), "1 categories")

  # Exemplar model with n_exemplars = 0 (default) vs n_exemplars > 0
  ex_df <- data.frame(F1 = seq(100, 500, length.out = 30), category = "A")
  ex_rep <- new_exemplar_category_representation_from_data(ex_df, cues = "F1", category = "category")
  ex_tpl <- new_category_representation_template(list(A = ex_rep))

  p_ex_default <- plot_categories(ex_tpl)
  expect_no_match(deparse1(p_ex_default$labels$subtitle), "sampling .* exemplars")

  p_ex_sampled <- plot_categories(ex_tpl, n_exemplars = 15L)
  expect_match(deparse1(p_ex_sampled$labels$subtitle), "sampling 15 exemplars")
})

test_that("3D interactive category plots work with plotly and various aes options", {
  skip_if_not_installed("plotly")

  # Parametric 3D model
  tpl3d <- example_mvg_category_representation_template(n_cues = 3)
  p3d_param <- plot_categories(tpl3d, interactive = TRUE)
  expect_s3_class(p3d_param, "plotly")

  # Parametric with contour aes
  p3d_cont <- plot_categories(tpl3d, interactive = TRUE, aes = "contour")
  expect_s3_class(p3d_cont, "plotly")

  # Parametric with fill-gradient aes
  p3d_grad <- plot_categories(tpl3d, interactive = TRUE, aes = "fill-gradient")
  expect_s3_class(p3d_grad, "plotly")

  # Exemplar 3D model with scatter
  ex3d <- example_exemplar_category_representation_template(n_cues = 3)
  p3d_ex <- plot_categories(ex3d, aes = "scatter", interactive = TRUE, n_exemplars = 50L)
  expect_s3_class(p3d_ex, "plotly")

  # Exemplar 3D model with isosurface fill
  p3d_ex_iso <- plot_categories(ex3d, aes = "fill-discrete", interactive = TRUE)
  expect_s3_class(p3d_ex_iso, "plotly")
})

test_that("2D interactive category plots render 3D density surfaces via plotly", {
  skip_if_not_installed("plotly")

  tpl2d <- example_mvg_category_representation_template(n_cues = 2)
  p2d_disc <- plot_categories(tpl2d, interactive = TRUE, aes = "fill-discrete")
  expect_s3_class(p2d_disc, "plotly")

  p2d_grad <- plot_categories(tpl2d, interactive = TRUE, aes = "fill-gradient")
  expect_s3_class(p2d_grad, "plotly")

  p2d_cont <- plot_categories(tpl2d, interactive = TRUE, aes = "contour")
  expect_s3_class(p2d_cont, "plotly")
})

test_that("2D exemplar models support fill-discrete aesthetic", {
  ex_df <- data.frame(
    F1 = rnorm(50, 300, 30),
    F2 = rnorm(50, 1000, 80),
    category = "A"
  )
  ex_rep <- new_exemplar_category_representation_from_data(
    ex_df,
    cues = c("F1", "F2"),
    category = "category"
  )
  ex_tpl <- new_category_representation_template(list(A = ex_rep))

  p_disc <- plot_categories(ex_tpl, aes = "fill-discrete")
  expect_s3_class(p_disc, "ggplot")
  layer_types <- vapply(p_disc$layers, function(l) class(l$geom)[1], character(1))
  expect_true("GeomPolygon" %in% layer_types)
})

test_that("3D sliced category plots accept levels and slices aliases and default to 5 slices", {
  tpl3d <- example_mvg_category_representation_template(n_cues = 3)
  cues3 <- get_cue_labels(tpl3d)[1:3]

  # Default slices (-2, -1, 0, 1, 2 sigma) and default levels (1:3 sigma)
  p_def <- plot_categories(tpl3d, cues = cues3, aes = c("fill-gradient", "contour"))
  expect_s3_class(p_def, "ggplot")
  layer_types <- vapply(p_def$layers, function(l) class(l$geom)[1], character(1))
  expect_true("GeomPath" %in% layer_types || "GeomText" %in% layer_types)

  # Custom slices
  p_sliced <- plot_categories(
    tpl3d,
    cues = cues3,
    slices = c(0, 10, 20),
    levels = c(0.68, 0.95),
    aes = c("fill-gradient", "contour")
  )
  expect_s3_class(p_sliced, "ggplot")
})

test_that("cue defaulting selects up to 3 cues and error triggers when > 3 cues", {
  # 4D MVG model
  rep4d_A <- new_mvg_category_representation(
    category_labels = "A",
    cue_labels = c("c1", "c2", "c3", "c4"),
    mu = c(0, 1, 2, 3),
    Sigma = diag(4)
  )
  rep4d_B <- new_mvg_category_representation(
    category_labels = "B",
    cue_labels = c("c1", "c2", "c3", "c4"),
    mu = c(1, 2, 3, 4),
    Sigma = diag(4)
  )
  tpl4d <- new_category_representation_template(list(A = rep4d_A, B = rep4d_B))
  m4d <- new_mvg_ideal_observer(
    category_template = tpl4d,
    category_prior = c(A = 0.5, B = 0.5)
  )

  # plot_categories defaults to first 3 cues (sliced plot)
  p_cat_def <- plot_categories(m4d)
  expect_s3_class(p_cat_def, "ggplot")

  # plot_categorization_functions defaults to first 3 cues (sliced plot)
  p_func_def <- plot_categorization_functions(m4d, categories = "A", resolution = 15)
  expect_s3_class(p_func_def, "ggplot")

  # Error when > 3 cues passed explicitly
  expect_error(
    plot_categories(m4d, cues = c("c1", "c2", "c3", "c4")),
    "Cannot plot more than 3"
  )
  expect_error(
    plot_categorization_functions(m4d, cues = c("c1", "c2", "c3", "c4")),
    "Cannot plot more than 3"
  )
})

test_that("2D interactive categorization plot returns plotly widget", {
  tpl2d <- example_mvg_category_representation_template(n_cues = 2)
  cats <- get_category_labels(tpl2d)
  cues <- get_cue_labels(tpl2d)
  m2d <- new_mvg_ideal_observer(
    category_template = tpl2d,
    category_prior = stats::setNames(rep(1 / length(cats), length(cats)), cats)
  )

  p_inter <- plot_categorization_functions(
    m2d,
    cues = cues,
    categories = cats[1L],
    interactive = TRUE,
    resolution = 15
  )
  expect_s3_class(p_inter, "plotly")
})

test_that("plot_sample and shorthands work on 1D and 2D inputs and fits", {
  set.seed(42)
  inp1d <- new_ideal_adaptor_stanfit_input(
    exposure = data.frame(
      category = factor(rep(c("A", "B"), each = 20)),
      group = factor(rep(c("g1", "g2"), 20)),
      cue1 = rnorm(40, mean = rep(c(0, 2), each = 20))
    ),
    test = data.frame(
      response = factor(rep(c("A", "B"), each = 10)),
      group = factor(rep(c("g1", "g2"), 10)),
      cue1 = rnorm(20, mean = 1)
    ),
    cues = "cue1",
    category = "category",
    response = "response",
    group = "group",
    stanmodel = "NIW_ideal_adaptor"
  )

  # 1D plot_sample, plot_exposure_sample, and plot_test_sample
  p_1d <- plot_sample(inp1d)
  expect_s3_class(p_1d, "ggplot")

  p_1d_exp <- plot_exposure_sample(inp1d)
  expect_s3_class(p_1d_exp, "ggplot")

  p_1d_tst <- plot_test_sample(inp1d)
  expect_s3_class(p_1d_tst, "ggplot")

  # 1D densities: gaussian vs kernel
  p_1d_gauss <- plot_exposure_sample(inp1d, densities = "gaussian")
  expect_s3_class(p_1d_gauss, "ggplot")

  p_1d_kern <- plot_exposure_sample(inp1d, densities = "kernel")
  expect_s3_class(p_1d_kern, "ggplot")

  # 2D plot_sample with various aesthetics
  inp2d <- new_ideal_adaptor_stanfit_input(
    exposure = data.frame(
      category = factor(rep(c("A", "B"), each = 20)),
      group = factor(rep(c("g1", "g2"), 20)),
      Condition = factor(rep(c("c1", "c2"), 20)),
      cue1 = rnorm(40, mean = rep(c(0, 2), each = 20)),
      cue2 = rnorm(40, mean = rep(c(0, 2), each = 20))
    ),
    test = data.frame(
      response = factor(rep(c("A", "B"), each = 10)),
      group = factor(rep(c("g1", "g2"), 10)),
      Condition = factor(rep(c("c1", "c2"), 10)),
      cue1 = rnorm(20, mean = 1),
      cue2 = rnorm(20, mean = 1)
    ),
    cues = c("cue1", "cue2"),
    category = "category",
    response = "response",
    group = "group",
    group_unique = "Condition",
    stanmodel = "NIW_ideal_adaptor"
  )

  p_2d_def <- plot_sample(inp2d)
  expect_s3_class(p_2d_def, "ggplot")

  p_2d_pts <- plot_sample(inp2d, aes = "points")
  expect_s3_class(p_2d_pts, "ggplot")

  p_2d_cont <- plot_sample(inp2d, aes = "contour")
  expect_s3_class(p_2d_cont, "ggplot")

  p_2d_grad <- plot_sample(inp2d, aes = "fill-gradient")
  expect_s3_class(p_2d_grad, "ggplot")

  p_2d_disc <- plot_sample(inp2d, aes = "fill-discrete")
  expect_s3_class(p_2d_disc, "ggplot")

  # aes list specification
  p_2d_list_aes <- plot_sample(
    inp2d,
    aes = list(exposure = c("fill-gradient", "contour"), test = "points")
  )
  expect_s3_class(p_2d_list_aes, "ggplot")

  # Subsampling with n.samples
  p_2d_sub <- plot_sample(inp2d, n.samples = 15)
  expect_s3_class(p_2d_sub, "ggplot")

  # Filtering by category and group
  p_2d_filtered <- plot_sample(
    inp2d,
    categories = get_category_labels(inp2d)[1],
    groups = get_group_labels(inp2d)[1]
  )
  expect_s3_class(p_2d_filtered, "ggplot")

  # Test points styling has distinct/black color
  p_test_pts <- plot_test_sample(inp2d, aes = "points")
  expect_s3_class(p_test_pts, "ggplot")

  # "all" is no longer valid; c("exposure", "test") is the default
  expect_error(plot_sample(inp2d, sample = "all"))
  p_both <- plot_sample(inp2d, sample = c("exposure", "test"))
  expect_s3_class(p_both, "ggplot")
  p_sub <- plot_sample(inp2d, n_samples = 35)
  expect_s3_class(p_sub, "ggplot")

  # 3D plot_sample and plot_exposure_sample
  inp3d <- new_ideal_adaptor_stanfit_input(
    exposure = data.frame(
      category = factor(rep(c("A", "B"), each = 25)),
      group = factor(rep("g1", 50)),
      cue1 = rnorm(50, 0, 1),
      cue2 = rnorm(50, 2, 1),
      cue3 = rnorm(50, -1, 1)
    ),
    test = data.frame(
      response = factor(rep(c("A", "B"), each = 15)),
      group = factor(rep("g1", 30)),
      cue1 = rnorm(30, 0, 1),
      cue2 = rnorm(30, 2, 1),
      cue3 = rnorm(30, -1, 1)
    ),
    cues = c("cue1", "cue2", "cue3"),
    category = "category",
    response = "response",
    group = "group",
    stanmodel = "NIW_ideal_adaptor"
  )

  p_3d <- plot_sample(inp3d, cues = c("cue1", "cue2", "cue3"))
  expect_s3_class(p_3d, "ggplot")

  p_3d_exp <- plot_exposure_sample(inp3d, cues = c("cue1", "cue2", "cue3"))
  expect_s3_class(p_3d_exp, "ggplot")

  p_3d_kern <- plot_exposure_sample(inp3d, cues = c("cue1", "cue2", "cue3"), densities = "kernel")
  expect_s3_class(p_3d_kern, "ggplot")

  # Error on non-existent cue
  expect_error(plot_sample(inp2d, cues = "NonExistentCue"))

  # Error on > 3 cues
  inp4d <- new_ideal_adaptor_stanfit_input(
    exposure = data.frame(
      category = factor(rep(c("A", "B"), each = 20)),
      group = factor(rep("g1", 40)),
      c1 = rnorm(40), c2 = rnorm(40), c3 = rnorm(40), c4 = rnorm(40)
    ),
    test = data.frame(
      response = factor(rep(c("A", "B"), each = 10)),
      group = factor(rep("g1", 20)),
      c1 = rnorm(20), c2 = rnorm(20), c3 = rnorm(20), c4 = rnorm(20)
    ),
    cues = c("c1", "c2", "c3", "c4"),
    category = "category",
    response = "response",
    group = "group",
    stanmodel = "NIW_ideal_adaptor"
  )
  expect_error(
    plot_sample(inp4d, cues = c("c1", "c2", "c3", "c4")),
    "Cannot plot more than 3 cue dimensions simultaneously"
  )
})

test_that("get_original_variable_names and get_data(original_names = TRUE) work", {
  inp <- new_ideal_adaptor_stanfit_input(
    exposure = data.frame(
      category = factor(c("A", "B")),
      group = factor(c("g1", "g2")),
      cue1 = c(1, 2)
    ),
    test = data.frame(
      response = factor(c("A", "B")),
      group = factor(c("g1", "g2")),
      cue1 = c(1.5, 2.5)
    ),
    cues = "cue1",
    category = "category",
    response = "response",
    group = "group",
    stanmodel = "NIW_ideal_adaptor"
  )

  # Default returns full named list
  orig_all <- get_original_variable_names(inp)
  expect_true(is.list(orig_all))
  expect_true(all(c("group", "group_unique", "category", "response_category", "cues") %in% names(orig_all)))

  # Querying single variable returns character vector
  orig_cat <- get_original_variable_names(inp, "category")
  expect_true(is.character(orig_cat))
  expect_equal(length(orig_cat), 1L)

  orig_cues <- get_original_variable_names(inp, "cues")
  expect_true(is.character(orig_cues))

  # Querying multiple variables returns subset list
  orig_sub <- get_original_variable_names(inp, c("group", "category"))
  expect_true(is.list(orig_sub))
  expect_equal(names(orig_sub), c("group", "category"))

  # get_data with original_names = FALSE (default) returns standard column names
  std_data <- get_data(inp, original_names = FALSE)
  expect_true("category" %in% names(std_data))
  expect_true("response_category" %in% names(std_data))
  expect_true("group" %in% names(std_data))
  expect_true("group_unique" %in% names(std_data))

  # get_data with original_names = TRUE returns user's original column names
  orig_data <- get_data(inp, original_names = TRUE)
  expect_true(orig_all$group %in% names(orig_data))
  expect_true(orig_all$category %in% names(orig_data))
  expect_true(orig_all$response_category %in% names(orig_data))

  # Same for get_exposure_data and get_test_data
  exp_std <- get_exposure_data(inp, original_names = FALSE)
  expect_true("category" %in% names(exp_std))
  exp_orig <- get_exposure_data(inp, original_names = TRUE)
  expect_true(orig_all$category %in% names(exp_orig))

  tst_std <- get_test_data(inp, original_names = FALSE)
  expect_true("response_category" %in% names(tst_std))
  tst_orig <- get_test_data(inp, original_names = TRUE)
  expect_true(orig_all$response_category %in% names(tst_orig))
})


