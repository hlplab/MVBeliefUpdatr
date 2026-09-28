test_that("as_tibble on MVG template and model yields exact legacy schema and classes", {
  mvg_model <- example_mvg_ideal_observer(n_cues = 2)
  mvg_tpl <- mvg_model@category_template

  # Template
  expect_warning(
    tpl_df <- tibble::as_tibble(mvg_tpl),
    "deprecated"
  )
  expect_s3_class(tpl_df, "tbl_df")
  expect_named(tpl_df, c("category", "mu", "Sigma"))
  expect_true(is.factor(tpl_df$category))
  expect_equal(levels(tpl_df$category), c("/b/", "/d/", "/g/", "/k/", "/p/", "/t/"))
  expect_equal(names(tpl_df$mu[[1]]), c("VOT", "f0_semitones"))
  expect_equal(dimnames(tpl_df$Sigma[[1]]), list(c("VOT", "f0_semitones"), c("VOT", "f0_semitones")))
  expect_true(suppressWarnings(is.MVG(tpl_df)))
  expect_true(suppressWarnings(is.MVBU_representation(tpl_df)))

  # Model
  expect_warning(
    model_df <- tibble::as_tibble(mvg_model),
    "deprecated"
  )
  expect_s3_class(model_df, "tbl_df")
  expect_named(model_df, c("category", "mu", "Sigma", "prior", "lapse_rate", "lapse_bias", "Sigma_noise"))
  expect_equal(model_df$prior, get_category_prior(mvg_model))
  expect_equal(model_df$lapse_rate, rep(get_lapse_rate(mvg_model), nrow(model_df)))
  expect_equal(model_df$lapse_bias, get_lapse_bias(mvg_model))
  expect_true(suppressWarnings(is.MVG_ideal_observer(model_df)))
  expect_true(suppressWarnings(is.MVBU_model(model_df)))
})

test_that("as_tibble on NIW template and model yields exact legacy schema and classes", {
  niw_model <- example_niw_ideal_adaptor(n_cues = 2)
  niw_tpl <- niw_model@category_template

  # Template
  expect_warning(
    tpl_df <- tibble::as_tibble(niw_tpl),
    "deprecated"
  )
  expect_s3_class(tpl_df, "tbl_df")
  expect_named(tpl_df, c("category", "m", "S", "kappa", "nu"))
  expect_true(is.factor(tpl_df$category))
  expect_equal(names(tpl_df$m[[1]]), c("VOT", "f0_semitones"))
  expect_equal(dimnames(tpl_df$S[[1]]), list(c("VOT", "f0_semitones"), c("VOT", "f0_semitones")))
  expect_true(suppressWarnings(is.NIW_belief(tpl_df)))
  expect_true(suppressWarnings(is.MVBU_representation(tpl_df)))

  # Model
  expect_warning(
    model_df <- tibble::as_tibble(niw_model),
    "deprecated"
  )
  expect_s3_class(model_df, "tbl_df")
  expect_named(model_df, c("category", "m", "S", "kappa", "nu", "prior", "lapse_rate", "lapse_bias", "Sigma_noise"))
  expect_true(suppressWarnings(is.NIW_ideal_adaptor(model_df)))
  expect_true(suppressWarnings(is.MVBU_model(model_df)))
})

test_that("as_tibble on Exemplar template and model yields exact legacy schema and classes", {
  ex_model <- example_exemplar_model(n_cues = 2)
  ex_tpl <- ex_model@category_template

  # Template
  expect_warning(
    tpl_df <- tibble::as_tibble(ex_tpl),
    "deprecated"
  )
  expect_s3_class(tpl_df, "tbl_df")
  expect_named(tpl_df, c("category", "exemplars", "sim_function"))
  expect_true(is.matrix(tpl_df$exemplars[[1]]))
  expect_true(is.function(tpl_df$sim_function[[1]]))
  expect_true(suppressWarnings(is.exemplars(tpl_df)))
  expect_true(suppressWarnings(is.MVBU_representation(tpl_df)))

  # Model
  expect_warning(
    model_df <- tibble::as_tibble(ex_model),
    "deprecated"
  )
  expect_s3_class(model_df, "tbl_df")
  expect_named(model_df, c("category", "exemplars", "sim_function", "prior", "lapse_rate", "lapse_bias", "Sigma_noise"))
  expect_true(suppressWarnings(is.exemplar_model(model_df)))
  expect_true(suppressWarnings(is.MVBU_model(model_df)))
})

test_that("as_tibble on UVG and NIX templates and models", {
  uvg_model <- example_uvg_ideal_observer()
  expect_warning(
    uvg_df <- tibble::as_tibble(uvg_model),
    "deprecated"
  )
  expect_s3_class(uvg_df, "tbl_df")
  expect_named(uvg_df, c("category", "mu", "sigma", "prior", "lapse_rate", "lapse_bias", "Sigma_noise"))

  nix_model <- example_nix_ideal_adaptor()
  expect_warning(
    nix_df <- tibble::as_tibble(nix_model),
    "deprecated"
  )
  expect_s3_class(nix_df, "tbl_df")
  expect_named(nix_df, c("category", "m", "S", "kappa", "nu", "prior", "lapse_rate", "lapse_bias", "Sigma_noise"))
})

test_that("as_tibble on single category representation works", {
  mvg_model <- example_mvg_ideal_observer(n_cues = 2)
  rep1 <- mvg_model@category_template@representations[[1]]

  expect_warning(
    df <- tibble::as_tibble(rep1),
    "deprecated"
  )
  expect_s3_class(df, "tbl_df")
  expect_equal(nrow(df), 1L)
  expect_named(df, c("category", "mu", "Sigma"))
})

test_that("as.data.frame on S7 model delegates through as_tibble", {
  mvg_model <- example_mvg_ideal_observer(n_cues = 2)
  expect_warning(
    df <- as.data.frame(mvg_model),
    "deprecated"
  )
  expect_true(is.data.frame(df))
  expect_true("prior" %in% names(df))
})

test_that("deprecated dplyr verbs forward seamlessly on S7 objects", {
  mvg_model <- example_mvg_ideal_observer(n_cues = 2)

  # mutate
  expect_warning(
    res_mutate <- dplyr::mutate(mvg_model, double_prior = prior * 2),
    "deprecated"
  )
  expect_s3_class(res_mutate, "tbl_df")
  expect_true("double_prior" %in% names(res_mutate))
  expect_equal(res_mutate$double_prior, res_mutate$prior * 2)

  # filter
  expect_warning(
    res_filter <- dplyr::filter(mvg_model, category == "/b/"),
    "deprecated"
  )
  expect_s3_class(res_filter, "tbl_df")
  expect_equal(nrow(res_filter), 1L)
  expect_equal(as.character(res_filter$category), "/b/")

  # select
  expect_warning(
    res_select <- dplyr::select(mvg_model, category, prior),
    "deprecated"
  )
  expect_named(res_select, c("category", "prior"))

  # pull
  expect_warning(
    priors <- dplyr::pull(mvg_model, prior),
    "deprecated"
  )
  expect_equal(priors, get_category_prior(mvg_model))

  # arrange
  expect_warning(
    res_arrange <- dplyr::arrange(mvg_model, dplyr::desc(category)),
    "deprecated"
  )
  expect_s3_class(res_arrange, "tbl_df")

  # transmute
  expect_warning(
    res_transmute <- dplyr::transmute(mvg_model, category, p2 = prior * 2),
    "deprecated"
  )
  expect_named(res_transmute, c("category", "p2"))

  # slice
  expect_warning(
    res_slice <- dplyr::slice(mvg_model, 1:2),
    "deprecated"
  )
  expect_equal(nrow(res_slice), 2L)

  # rename
  expect_warning(
    res_rename <- dplyr::rename(mvg_model, cat = category),
    "deprecated"
  )
  expect_true("cat" %in% names(res_rename))

  # relocate
  expect_warning(
    res_reloc <- dplyr::relocate(mvg_model, prior, .before = category),
    "deprecated"
  )
  expect_equal(names(res_reloc)[1], "prior")

  # Pipe chaining: warning is only issued once at the first conversion
  expect_warning(
    piped <- mvg_model %>%
      dplyr::filter(prior > 0) %>%
      dplyr::mutate(test_col = 1),
    "deprecated"
  )
  expect_true("test_col" %in% names(piped))
})
