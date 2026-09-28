test_that("deprecated get_constructor functions raise lifecycle warnings", {
  skip_if_not_installed("rstan")
  skip_if_not_installed("tidybayes")

  model_path <- testthat::test_path(
    "models",
    "minimal-nix_ideal_adaptor-NIX.rds"
  )
  fit <- example_ideal_adaptor_stanfit(
    model_family = "NIX",
    n_cues = 1L,
    file = model_path,
    file_refit = "on_change"
  )

  recovered_stanfit <- tidybayes::recover_types(get_stanfit(fit))
  recovered_fit <- set_stanfit(fit, recovered_stanfit)

  expect_warning(
    ctor_group <- get_group_constructor(recovered_fit),
    class = "lifecycle_warning_deprecated"
  )
  expect_true(is.function(ctor_group))

  expect_warning(
    ctor_cat <- get_category_constructor(recovered_fit),
    class = "lifecycle_warning_deprecated"
  )
  expect_true(is.function(ctor_cat))

  expect_warning(
    ctor_all <- get_constructor(recovered_fit),
    class = "lifecycle_warning_deprecated"
  )
  expect_type(ctor_all, "list")

  expect_warning(
    ctor_cue <- get_cue_constructor(recovered_fit),
    class = "lifecycle_warning_deprecated"
  )
  expect_true(is.function(ctor_cue))

  expect_warning(
    ctor_cue2 <- get_cue2_constructor(recovered_fit),
    class = "lifecycle_warning_deprecated"
  )
})
