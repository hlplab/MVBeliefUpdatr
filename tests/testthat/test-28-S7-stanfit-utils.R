test_that(".check_stanfit_file normalizes file extensions correctly", {
  expect_equal(MVBeliefUpdatr:::.check_stanfit_file("my_model"), "my_model.rds")
  expect_equal(MVBeliefUpdatr:::.check_stanfit_file("my_model.rds"), "my_model.rds")
  expect_equal(MVBeliefUpdatr:::.check_stanfit_file("path/to/fit.RDS"), "path/to/fit.RDS")
  expect_error(MVBeliefUpdatr:::.check_stanfit_file(NA_character_))
  expect_error(MVBeliefUpdatr:::.check_stanfit_file(c("a", "b")))
})

test_that(".file_refit_options returns expected choices", {
  expect_equal(
    MVBeliefUpdatr:::.file_refit_options(),
    c("never", "always", "on_change")
  )
})

test_that("write_stanfit and read_stanfit roundtrip an MVBU_Stanfit object", {
  fit <- get_example_stanfit(1, file_refit = "never")
  tmp_dir <- withr::local_tempdir()
  tmp_path <- file.path(tmp_dir, "test_fit") # test auto-appending .rds

  # write_stanfit
  res <- write_stanfit(fit, tmp_path)
  expected_file <- paste0(tmp_path, ".rds")
  expect_true(file.exists(expected_file))
  expect_equal(res@file, expected_file)

  # read_stanfit
  loaded <- read_stanfit(tmp_path)
  expect_true(S7::S7_inherits(loaded, MVBU_Stanfit))
  expect_true(S7::S7_inherits(loaded, IdealAdaptorStanfit))
  expect_equal(loaded@file, expected_file)
  expect_equal(loaded@version, fit@version)

  # read_stanfit returns NULL on non-existent file
  expect_null(read_stanfit(file.path(tmp_dir, "nonexistent")))

  # read_stanfit returns NULL on corrupted file
  bad_file <- file.path(tmp_dir, "corrupt.rds")
  writeLines("not an rds", bad_file)
  expect_null(read_stanfit(bad_file))

  # write_stanfit fails on non-MVBU_Stanfit object
  expect_error(write_stanfit(list(a = 1), file.path(tmp_dir, "fail.rds")))
})

test_that(".stanfit_needs_refit behaves correctly", {
  fit <- get_example_stanfit(1, file_refit = "never")

  # Current version and unchanged inputs -> no refit
  expect_false(
    MVBeliefUpdatr:::.stanfit_needs_refit(
      fit,
      current_version = fit@version,
      silent = TRUE
    )
  )

  # Changed package version -> refit needed
  diff_version <- fit@version
  diff_version$MVBeliefUpdatr <- "99.99.99"
  expect_true(
    MVBeliefUpdatr:::.stanfit_needs_refit(
      fit,
      current_version = diff_version,
      silent = TRUE
    )
  )

  # Same staninput -> no refit
  cached_input <- get_staninput(fit)
  expect_false(
    MVBeliefUpdatr:::.stanfit_needs_refit(
      fit,
      current_version = fit@version,
      staninput = cached_input,
      silent = TRUE
    )
  )

  # Modified staninput -> refit needed
  modified_input <- cached_input
  if (S7::S7_inherits(modified_input, IdealAdaptorStaninput)) {
    modified_input@values[[1]] <- modified_input@values[[1]] + 1
  } else if (is.list(modified_input) && length(modified_input) > 0) {
    modified_input[[1]] <- modified_input[[1]] + 1
  }
  expect_true(
    MVBeliefUpdatr:::.stanfit_needs_refit(
      fit,
      current_version = fit@version,
      staninput = modified_input,
      silent = TRUE
    )
  )

  # Data with identical factor levels -> no refit
  if (!is.null(fit@data) && nrow(fit@data) > 0) {
    expect_false(
      MVBeliefUpdatr:::.stanfit_needs_refit(
        fit,
        current_version = fit@version,
        data = fit@data,
        silent = TRUE
      )
    )

    # Data with changed factor levels -> refit needed
    mod_data <- fit@data
    cat_col <- which(vapply(mod_data, is.factor, logical(1)))[1]
    if (!is.na(cat_col)) {
      levels(mod_data[[cat_col]]) <- rev(levels(mod_data[[cat_col]]))
      expect_true(
        MVBeliefUpdatr:::.stanfit_needs_refit(
          fit,
          current_version = fit@version,
          data = mod_data,
          silent = TRUE
        )
      )
    }
  }
})
