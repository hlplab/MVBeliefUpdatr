# Tests for deprecated is.ideal_adaptor_stanfit() type check function.
# These tests verify the deprecated function still works correctly while
# emitting lifecycle deprecation warnings.

skip_if_not(file.exists(example_stanfit_path(1, stanmodel = "NIW_ideal_adaptor")), "cached example stanfit not generated")
fit <- get_example_stanfit(1, stanmodel = "NIW_ideal_adaptor", file_refit = "never")

test_that("Deprecated is.ideal_adaptor_stanfit emits warning and returns correct values", {
  suppressWarnings({
    expect_false(is.ideal_adaptor_stanfit(NULL))
    expect_false(is.ideal_adaptor_stanfit(NA))
    expect_false(is.ideal_adaptor_stanfit(1))
    expect_false(is.ideal_adaptor_stanfit("1"))
    expect_false(is.ideal_adaptor_stanfit(TRUE))
    expect_false(is.ideal_adaptor_stanfit(list(1)))
    expect_false(is.ideal_adaptor_stanfit(example_exemplar_model(n_cues = 1)))
    expect_false(is.ideal_adaptor_stanfit(example_mvg_ideal_observer(n_cues = 1)))
    expect_false(is.ideal_adaptor_stanfit(example_niw_ideal_adaptor(n_cues = 1)))
    expect_true(is.ideal_adaptor_stanfit(fit))
  })

  # Verify the deprecation warning is emitted
  expect_warning(is.ideal_adaptor_stanfit(fit), "deprecated")
})
