test_that("update_template performs an S7 NIW batch update", {
  model <- example_niw_ideal_adaptor(n_cues = 2)
  categories <- get_category_labels(model)[1:2]
  exposure <- tibble::tibble(
    category = rep(categories, each = 2),
    VOT = c(10, 20, 40, 50),
    f0_semitones = c(-1, 0, 1, 2)
  )

  updated <- update_template(model, exposure)
  before <- model@category_template@representations
  after <- updated@category_template@representations

  expect_true(S7::S7_inherits(updated, NIW_IdealAdaptor))
  expect_equal(after[[1]]@kappa, before[[1]]@kappa + 2)
  expect_equal(after[[1]]@nu, before[[1]]@nu + 2)
  expect_false(identical(after[[1]]@m, before[[1]]@m))
  expect_false(identical(after[[1]]@S, before[[1]]@S))
  expect_equal(after[[2]]@kappa, before[[2]]@kappa + 2)
  expect_equal(after[[2]]@nu, before[[2]]@nu + 2)
})

test_that("update_category_representation updates an NIW representation", {
  representation <- example_niw_category_representation(n_cues = 2, nu = 30)
  updated <- update_category_representation(
    representation, x_N = 2, x_mean = c(15, .5),
    x_SS = matrix(c(50, 5, 5, .5), nrow = 2)
  )
  expect_true(S7::S7_inherits(updated, NIW_CategoryRepresentation))
  expect_equal(updated@kappa, representation@kappa + 2)
  expect_equal(updated@nu, representation@nu + 2)
  expect_false(identical(updated@m, representation@m))
  expect_false(identical(updated@S, representation@S))
})

test_that("update_template preserves categories without exposure", {
  model <- example_niw_ideal_adaptor(n_cues = 1)
  category <- get_category_labels(model)[1]
  exposure <- tibble::tibble(category = category, VOT = c(10, 20))
  updated <- update_template(model, exposure)
  expect_identical(
    updated@category_template@representations[[2]]@m,
    model@category_template@representations[[2]]@m
  )
  expect_identical(
    updated@category_template@representations[[2]]@S,
    model@category_template@representations[[2]]@S
  )
})

test_that("update_template rejects incompatible S7 update data", {
  model <- example_niw_ideal_adaptor(n_cues = 1)
  expect_error(update_template(model, data.frame(category = "A", wrong_cue = 1)))
  expect_error(update_template(model, data.frame(wrong_category = "A", VOT = 1)))
})

test_that("update_template supports incremental history and update methods", {
  model <- example_niw_ideal_adaptor(n_cues = 1)
  observations <- tibble::tibble(
    category = get_category_labels(model)[1],
    VOT = 10
  )

  history <- update_template(
    model,
    observations,
    updating = "incremental",
    keep_history = TRUE,
    lapse_treatment = "no_lapses",
    noise_treatment = "no_noise",
    update_method = "label-certain"
  )
  expect_type(history, "list")
  expect_length(history, 2)
  expect_true(all(vapply(history, function(x) S7::S7_inherits(x, NIW_IdealAdaptor), logical(1))))

  unlabeled <- observations[0, c("VOT"), drop = FALSE]
  expect_true(S7::S7_inherits(
    update_template(model, unlabeled, updating = "incremental", update_method = "nolabel-uniform"),
    NIW_IdealAdaptor
  ))
  expect_error(update_template(model, observations, updating = "batch", update_method = "nolabel-uniform"))
  expect_error(update_template(model, observations, update_method = "not-a-method"))
})
