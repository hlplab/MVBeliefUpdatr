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
  fit <- get_example_stanfit(1)
  tmp_dir <- tempfile("mvbu_test_")
  dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
  withr::defer(unlink(tmp_dir, recursive = TRUE))
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
  fit <- get_example_stanfit(1)

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
  fit_data <- get_data(fit, space = "original")
  if (!is.null(fit_data) && nrow(fit_data) > 0) {
    expect_false(
      MVBeliefUpdatr:::.stanfit_needs_refit(
        fit,
        current_version = fit@version,
        data = fit_data,
        silent = TRUE
      )
    )

    # Data with changed factor levels -> refit needed
    mod_data <- fit_data
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

  # Missing metadata$original_variable_names -> refit needed
  fit_no_meta <- fit
  fit_no_meta@metadata$original_variable_names <- NULL
  expect_true(
    MVBeliefUpdatr:::.stanfit_needs_refit(
      fit_no_meta,
      current_version = fit@version,
      silent = TRUE
    )
  )

  # Missing standard column from @data -> refit needed
  for (col in c(
    "category", "response_category", "group", "group_unique", "Phase"
  )) {
    fit_missing_col <- fit
    fit_missing_col@data$original[[col]] <- NULL
    expect_true(
      MVBeliefUpdatr:::.stanfit_needs_refit(
        fit_missing_col,
        current_version = fit@version,
        silent = TRUE
      )
    )
  }
})

test_that(".validate_threads and .use_threading work as expected", {
  expect_null(MVBeliefUpdatr:::.validate_threads(NULL))
  expect_equal(MVBeliefUpdatr:::.validate_threads(1), list(threads = 1L))
  expect_equal(MVBeliefUpdatr:::.validate_threads(4L), list(threads = 4L))
  expect_equal(MVBeliefUpdatr:::.validate_threads(list(threads = 2)), list(threads = 2L))
  expect_error(MVBeliefUpdatr:::.validate_threads(0))
  expect_error(MVBeliefUpdatr:::.validate_threads(-1))
  expect_error(MVBeliefUpdatr:::.validate_threads("abc"))

  expect_false(MVBeliefUpdatr:::.use_threading(NULL))
  expect_false(MVBeliefUpdatr:::.use_threading(1))
  expect_true(MVBeliefUpdatr:::.use_threading(2))
  expect_true(MVBeliefUpdatr:::.use_threading(list(threads = 4)))
})

test_that(".normalize_stancode strips comments and whitespace", {
  raw_code <- "
    // single line comment
    parameters {
      /* multi
         line
         comment */
      real x; // end line
    }
  "
  norm <- MVBeliefUpdatr:::.normalize_stancode(raw_code)
  expect_false(grepl("//", norm))
  expect_false(grepl("/\\*", norm))
  expect_equal(norm, "parameters { real x; }")
})

test_that(".get_stan_file and .get_stan_include_paths locate Stan files", {
  f_niw <- MVBeliefUpdatr:::.get_stan_file("NIW_ideal_adaptor")
  expect_true(file.exists(f_niw))
  expect_true(grepl("NIW_ideal_adaptor\\.stan$", f_niw))

  inc_paths <- MVBeliefUpdatr:::.get_stan_include_paths()
  expect_true(length(inc_paths) > 0)
  expect_true(any(dir.exists(inc_paths)))

  expect_error(MVBeliefUpdatr:::.get_stan_file("nonexistent_stan_model_12345"))
})

test_that(".compile_model dispatches to cmdstanr when installed", {
  skip_if_not_installed("cmdstanr")
  skip_if(
    is.null(tryCatch(cmdstanr::cmdstan_path(), error = function(e) NULL)),
    "CmdStan is not installed"
  )

  stan_file <- MVBeliefUpdatr:::.get_stan_file("NIX_ideal_adaptor")
  # Use compile = FALSE to verify syntax check and backend dispatch without full C++ compilation
  mod <- MVBeliefUpdatr:::.compile_model(
    model = stan_file,
    backend = "cmdstanr",
    compile = FALSE,
    silent = 2
  )
  expect_s3_class(mod, "CmdStanModel")
  expect_true(mod$check_syntax())
})

test_that(".extract_fitted_parameter_names_from_stanfit extracts canonical names", {
  # Mock object that responds to names()
  names_stanfit <- c(
    "m_0[1,1]", "m_n[1,1]", "S_0[1,1,1]", "S_n[1,1,1]",
    "weights[1]", "weights[2]", "lapse_rate", "lapse_bias[1]",
    "prior_p[1]", "prior_p[2]", "Sigma_noise[1,1]",
    "m_0_tau[1]", "m_0_param[1]", "weights_param[1]",
    "lp__", "log_lik[1]"
  )
  extracted <- MVBeliefUpdatr:::.extract_fitted_parameter_names_from_stanfit(
    names_stanfit
  )
  expect_true("m" %in% extracted)
  expect_true("S" %in% extracted)
  expect_true("weights" %in% extracted)
  expect_true("lapse_rate" %in% extracted)
  expect_true("lapse_bias" %in% extracted)
  expect_true("prior_p" %in% extracted)
  expect_true("Sigma_noise" %in% extracted)
  # Internal / Stan-specific names should not be extracted
  expect_false("lp__" %in% extracted)
  expect_false("log_lik" %in% extracted)
  expect_false("m_0_tau" %in% extracted)
  expect_false("cue_weight" %in% extracted)
  expect_false("p_category" %in% extracted)
})

test_that(".format_fixed_parameters_table formats parameters into standard tibble", {
  # Empty input
  empty_tbl <- MVBeliefUpdatr:::.format_fixed_parameters_table(list())
  expect_s3_class(empty_tbl, "tbl_df")
  expect_equal(nrow(empty_tbl), 0L)

  # NULL input
  expect_equal(
    nrow(MVBeliefUpdatr:::.format_fixed_parameters_table(NULL)),
    0L
  )

  # weights and lapse_rate
  fixed_pars <- list(
    weights = c(c1 = 0.7, c2 = 0.3),
    lapse_rate = 0.05,
    prior_p = c(catA = 0.6, catB = 0.4)
  )
  tbl <- MVBeliefUpdatr:::.format_fixed_parameters_table(fixed_pars)
  expect_s3_class(tbl, "tbl_df")
  expect_equal(nrow(tbl), 5L)
  expected_cols <- c(
    "Parameter", "Dist.", "Group", "Category", "Cue", "Value", "Status"
  )
  expect_true(all(expected_cols %in% names(tbl)))
  expect_equal(tbl$Status, rep("Fixed by user", 5L))

  w_rows <- tbl[tbl$Parameter == "weights", ]
  expect_equal(nrow(w_rows), 2L)
  expect_equal(w_rows$Cue, c("c1", "c2"))
  expect_equal(w_rows$Value, c(0.7, 0.3))

  lr_row <- tbl[tbl$Parameter == "lapse_rate", ]
  expect_equal(nrow(lr_row), 1L)
  expect_equal(lr_row$Value, 0.05)

  p_rows <- tbl[tbl$Parameter == "prior_p", ]
  expect_equal(nrow(p_rows), 2L)
  expect_equal(p_rows$Category, c("catA", "catB"))
  expect_equal(p_rows$Value, c(0.6, 0.4))
})



