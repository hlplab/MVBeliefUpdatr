# =============================================================================
# Unit Tests for MVBU_ModelList, MVBU_StanfitPosterior, & Multi-Model Plots
# =============================================================================

test_that("MVBU_ModelList creates, validates, and subsets correctly", {
  m1 <- example_nix_ideal_adaptor()
  m2 <- update_template(
    m1,
    observations = setNames(
      data.frame(1.0, "A"),
      c(get_cue_labels(m1)[1], "category")
    )
  )

  mod_list <- as_model_list(list(Prior = m1, Updated = m2))

  expect_true(S7::S7_inherits(mod_list, MVBU_ModelList))
  expect_equal(length(mod_list@models), 2L)
  expect_equal(mod_list@model_labels, c("Prior", "Updated"))

  # Subsetting tests
  sub_list <- mod_list[1]
  expect_true(S7::S7_inherits(sub_list, MVBU_ModelList))
  expect_equal(length(sub_list@models), 1L)
  expect_equal(sub_list@model_labels, "Prior")

  single_mod <- mod_list[[2]]
  expect_true(S7::S7_inherits(single_mod, NIX_IdealAdaptor))
})

test_that(
  "MVBU_StanfitPosterior creates lightweight representation from Stanfit",
  {
    sf_obj <- get_example_stanfit(
      1,
      stanmodel = "NIW_ideal_adaptor"
    )
    post_obj <- as_MVBU_stanfit_posterior(sf_obj)

    expect_true(S7::S7_inherits(post_obj, MVBU_StanfitPosterior))
    expect_equal(get_cue_labels(post_obj), get_cue_labels(sf_obj))
    expect_equal(get_category_labels(post_obj), get_category_labels(sf_obj))
    expect_equal(
      get_group_labels(post_obj),
      get_group_labels(sf_obj, include_prior = FALSE)
    )
    expect_equal(get_model_family(post_obj), get_model_family(sf_obj))
    expect_true(is.list(post_obj@draws))
    expect_true(is.list(get_metadata(post_obj)))
    expect_identical(
      post_obj@metadata$label_information$cue,
      get_cue_labels(sf_obj)
    )
    expect_identical(
      post_obj@metadata$label_information$category,
      get_category_labels(sf_obj)
    )

    suff <- get_sufficient_category_statistics(post_obj)
    expect_true(is.data.frame(suff))
    expect_false(inherits(suff, "tbl_df"))

    out_print <- utils::capture.output(print(post_obj))
    expect_true(any(grepl("Lightweight Stanfit Posterior", out_print)))
  }
)

test_that("reconstruct_update_history creates MVBU_ModelList from Stanfit", {
  sf_obj <- get_example_stanfit(
    1,
    stanmodel = "NIW_ideal_adaptor"
  )

  # uncertainty_treatment = "discard" (point estimate sequence)
  hist_fast <- reconstruct_update_history(
    sf_obj,
    uncertainty_treatment = "discard",
    step_size = 10L
  )
  expect_true(S7::S7_inherits(hist_fast, MVBU_ModelList))
  expect_true(length(hist_fast@models) >= 1L)

  # uncertainty_treatment = "marginalize" (MCMC posterior representation)
  hist_full <- reconstruct_update_history(
    sf_obj,
    uncertainty_treatment = "marginalize"
  )
  expect_true(S7::S7_inherits(hist_full, MVBU_ModelList))
  expect_true(S7::S7_inherits(hist_full[[1]], MVBU_StanfitPosterior))
})

test_that(
  "plot_categories and plot_categorization_functions work on MVBU_ModelList",
  {
    m1 <- example_nix_ideal_adaptor()
    m2 <- update_template(
      m1,
      observations = setNames(
        data.frame(1.0, "A"),
        c(get_cue_labels(m1)[1], "category")
      )
    )
    mod_list <- as_model_list(list(Prior = m1, Updated = m2))

    # Faceted ggplot outputs
    p_cat <- plot_categories(mod_list, mode = "facet", interactive = FALSE)
    expect_s3_class(p_cat, "ggplot")

    p_dec <- plot_categorization_functions(
      mod_list,
      mode = "facet",
      interactive = FALSE
    )
    expect_s3_class(p_dec, "ggplot")
  }
)

test_that(
  "plot_categories on MVBU_ModelList supports 2D single contour & 3D surface",
  {
    skip_if_not_installed("plotly")
    set.seed(42)
    m1_2d <- example_niw_ideal_adaptor(n_cues = 2)
    m2_2d <- update_template(m1_2d, data.frame(
      VOT = rnorm(15, mean = 25, sd = 5),
      f0_semitones = rnorm(15, mean = 3, sd = 1),
      category = "/b/"
    ))
    mod_list_2d <- as_model_list(list(Prior = m1_2d, Step1 = m2_2d))

    # 2D Plotly animation with single contour level 0.95
    p_2d <- plot_categories(
      mod_list_2d,
      mode = "animate",
      interactive = TRUE,
      levels = 0.95
    )
    expect_s3_class(p_2d, "plotly")
    pb_2d <- plotly::plotly_build(p_2d)
    expect_equal(length(pb_2d$x$frames), 2L)
    # 6 categories, 1 contour per category = 6 traces per frame
    expect_equal(length(pb_2d$x$frames[[1]]$data), 6L)

    # 3D Plotly animation with surface mesh for specific categories /b/ and /p/
    m1_3d <- example_niw_ideal_adaptor(n_cues = 3)
    m2_3d <- update_template(m1_3d, data.frame(
      VOT = rnorm(15, mean = 25, sd = 5),
      f0_semitones = rnorm(15, mean = 3, sd = 1),
      vowel_duration = rnorm(15, mean = 120, sd = 10),
      category = "/b/"
    ))
    mod_list_3d <- as_model_list(list(Prior = m1_3d, Step1 = m2_3d))

    p_3d <- plot_categories(
      mod_list_3d,
      mode = "animate",
      interactive = TRUE,
      categories = c("/b/", "/p/")
    )
    expect_s3_class(p_3d, "plotly")
    pb_3d <- plotly::plotly_build(p_3d)
    expect_equal(length(pb_3d$x$frames), 2L)
    expect_equal(length(pb_3d$x$frames[[1]]$data), 4L)
    types <- sapply(pb_3d$x$frames[[1]]$data, function(d) d$type)
    expect_equal(types, c("scatter3d", "mesh3d", "scatter3d", "mesh3d"))
    expect_true(length(pb_3d$x$layout$updatemenus) >= 1L)
    expect_equal(length(pb_3d$x$layout$updatemenus[[1]]$buttons), 1L)
    expect_equal(pb_3d$x$layout$updatemenus[[1]]$buttons[[1]]$label, "Play")
  }
)

test_that(
  "plot_categorization_functions on MVBU_ModelList supports 2D animated surface",
  {
    skip_if_not_installed("plotly")
    set.seed(42)
    m1_2d <- example_niw_ideal_adaptor(n_cues = 2)
    m2_2d <- update_template(m1_2d, data.frame(
      VOT = rnorm(15, mean = 25, sd = 5),
      f0_semitones = rnorm(15, mean = 3, sd = 1),
      category = "/b/"
    ))
    mod_list_2d <- as_model_list(list(Prior = m1_2d, Step1 = m2_2d))

    p_dec_anim <- plot_categorization_functions(
      mod_list_2d,
      mode = "animate",
      interactive = TRUE,
      categories = c("/b/", "/p/")
    )
    expect_s3_class(p_dec_anim, "plotly")
    pb_dec_anim <- plotly::plotly_build(p_dec_anim)
    expect_equal(length(pb_dec_anim$x$frames), 2L)
    expect_true(length(pb_dec_anim$x$layout$updatemenus) >= 1L)
    expect_equal(length(pb_dec_anim$x$layout$updatemenus[[1]]$buttons), 1L)
    expect_equal(pb_dec_anim$x$layout$updatemenus[[1]]$buttons[[1]]$label, "Play")
  }
)

