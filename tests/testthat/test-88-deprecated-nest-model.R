test_that("deprecated nest_cue_information_in_model and unnest_cue_information_in_model emit warnings", {
  # Create a test tibble with cue and cue2
  df_unnested <- tibble::tibble(
    group = rep("prior", 8),
    category = rep(c("/b/", "/d/"), each = 4),
    cue = rep(c("c1", "c1", "c2", "c2"), times = 2),
    cue2 = rep(c("c1", "c2", "c1", "c2"), times = 2),
    m = c(10, 10, 20, 20, 30, 30, 40, 40),
    S = c(1, 0, 0, 1, 2, 0, 0, 2)
  )

  expect_warning(
    df_nested <- nest_cue_information_in_model(df_unnested),
    class = "lifecycle_warning_deprecated"
  )
  expect_true("m" %in% names(df_nested))
  expect_true("S" %in% names(df_nested))
  expect_false("cue" %in% names(df_nested))

  expect_warning(
    df_restored <- unnest_cue_information_in_model(df_nested),
    class = "lifecycle_warning_deprecated"
  )
  expect_true("cue" %in% names(df_restored))
  expect_true("cue2" %in% names(df_restored))
})

test_that("deprecated make_named_vector and make_named_square_matrix emit warnings", {
  expect_warning(
    v <- make_named_vector(c(1, 2), c("a", "b")),
    class = "lifecycle_warning_deprecated"
  )
  expect_equal(names(v), c("a", "b"))

  expect_warning(
    mat <- make_named_square_matrix(c(1, 0, 0, 1), c("a", "b")),
    class = "lifecycle_warning_deprecated"
  )
  expect_equal(dim(mat), c(2, 2))
  expect_equal(rownames(mat), c("a", "b"))
})
