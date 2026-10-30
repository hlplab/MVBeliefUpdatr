test_that(
  "fit_ideal_adaptor fits and loads minimal examples across model families",
  {
    skip_if_not_installed("rstan")

    for (spec in list(
      list(
        family = "NIX",
        n_cues = 1L,
        file = "minimal-nix_ideal_adaptor-NIX.rds"
      ),
      list(
        family = "MNIX",
        n_cues = 2L,
        file = "minimal-mnix_ideal_adaptor-MNIX-2cue.rds"
      ),
      list(
        family = "NIW",
        n_cues = 2L,
        file = "minimal-niw_ideal_adaptor-NIW-2cue.rds"
      )
    )) {
      model_path <- testthat::test_path("models", spec$file)

      fit <- example_ideal_adaptor_stanfit(
        model_family = spec$family,
        n_cues = spec$n_cues,
        seed = 42L,
        file = model_path,
        file_refit = "on_change",
        chains = 1,
        iter = 500,
        warmup = 250,
        refresh = 0,
        control = list(adapt_delta = 0.99, max_treedepth = 12)
      )

      expect_true(
        file.exists(model_path),
        info = paste("model file was not written for", spec$family)
      )

      reloaded_fit <- read_stanfit(model_path)
      expect_true(
        S7::S7_inherits(reloaded_fit, IdealAdaptorStanfit),
        info = paste(
          "reloaded model is not an ideal adaptor fit for",
          spec$family
        )
      )

      expect_true(
        S7::S7_inherits(fit, IdealAdaptorStanfit),
        info = paste("fit object is not an ideal adaptor fit for", spec$family)
      )
      expect_s4_class(get_stanfit(fit), "stanfit")
      expect_true(
        .contains_draws(get_stanfit(fit)),
        info = paste("fit does not contain posterior draws for", spec$family)
      )
    }
  }
)


test_that(
  paste(
    "fit_ideal_adaptor loads an existing model from file when",
    "file_refit is never"
  ),
  {
    skip_if_not_installed("rstan")

    model_path <- testthat::test_path(
      "models",
      "minimal-nix_ideal_adaptor-NIX.rds"
    )
    fit <- read_stanfit(model_path)
    tmp_file <- tempfile(fileext = ".rds")
    write_stanfit(fit, tmp_file)

    input <- example_ideal_adaptor_stanfit_input(
      model_family = "NIX",
      n_cues = 1L,
      seed = 42L,
      control = control_staninput(transform_type = "identity")
    )

    reloaded_fit <- fit_ideal_adaptor(
      stanfit_input = input,
      file = tmp_file,
      file_refit = "never",
      chains = 0,
      iter = 0,
      warmup = 0,
      refresh = 0,
      stanmodel = "NIX_ideal_adaptor"
    )

  expect_true(S7::S7_inherits(reloaded_fit, IdealAdaptorStanfit))
  # stanfit@.MISC holds rstan's C++ module pointer, which is rebuilt on
  # every deserialization, so compare recovered content rather than whole obj.
  expect_equal(reloaded_fit@stanfit@sim$samples, fit@stanfit@sim$samples)
  expect_identical(reloaded_fit@stanfit@model_name, fit@stanfit@model_name)
  expect_identical(reloaded_fit@stanfit@model_pars, fit@stanfit@model_pars)
})

test_that(
  "recover_types works on stanfit objects in S7 IdealAdaptorStanfit objects",
  {
  skip_if_not_installed("rstan")
  skip_if_not_installed("tidybayes")

  model_path <- testthat::test_path(
    "models",
    "minimal-nix_ideal_adaptor-NIX.rds"
  )
  fit <- read_stanfit(model_path)

  recovered_stanfit <- tidybayes::recover_types(get_stanfit(fit))
  recovered_fit <- set_stanfit(fit, recovered_stanfit)

  expect_true(S7::S7_inherits(recovered_fit, IdealAdaptorStanfit))
  expect_s4_class(get_stanfit(recovered_fit), "stanfit")
  expect_warning(
    expect_true(is.function(get_constructor(recovered_fit, "group"))),
    class = "lifecycle_warning_deprecated"
  )
  expect_warning(
    val1 <- get_staninput_variable_levels(recovered_fit, "group"),
    class = "lifecycle_warning_deprecated"
  )
  expect_warning(
    val2 <- get_group_levels(recovered_fit),
    class = "lifecycle_warning_deprecated"
  )
  expect_equal(val1, val2)
})

test_that(
  "IdealAdaptorStaninput rejects empty family, preventing invalid inference",
  {
    skip_if_not_installed("rstan")
    raw_input <- example_ideal_adaptor_stanfit_input("NIW", n_cues = 2L)
    staninput <- MVBeliefUpdatr:::.new_staninput(
      raw_input,
      stanmodel = "NIW_ideal_adaptor"
    )

    expect_error(
      staninput@family <- "",
      regexp = paste0(
        "`family` must be a non-empty character string specifying a ",
        "registered Stan model family."
      )
    )
  }
)

test_that(
  "fit_ideal_adaptor throws an error on unrecognized model name",
  {
    raw_input <- example_ideal_adaptor_stanfit_input("MNIX", n_cues = 2L)
    expect_error(
      fit_ideal_adaptor(
        raw_input,
        stanmodel = "unknown_model_xyz"
      )
    )
  }
)

test_that("fit_ideal_adaptor supports backend = 'cmdstanr'", {
  skip_if_not_installed("cmdstanr")
  skip_if(
    is.null(tryCatch(cmdstanr::cmdstan_path(), error = function(e) NULL)),
    "CmdStan is not installed"
  )

  input <- example_ideal_adaptor_stanfit_input(
    model_family = "NIX",
    n_cues = 1L,
    seed = 42L,
    control = control_staninput(transform_type = "identity")
  )

  fit_empty <- fit_ideal_adaptor(
    stanfit_input = input,
    backend = "cmdstanr",
    chains = 0,
    iter = 0,
    warmup = 0,
    refresh = 0
  )
  expect_true(S7::S7_inherits(fit_empty, IdealAdaptorStanfit))
  expect_equal(fit_empty@backend, "cmdstanr")

  # Minimal sampling test to verify .read_cmdstan_csv_as_stanfit roundtrip
  fit_sample <- fit_ideal_adaptor(
    stanfit_input = input,
    backend = "cmdstanr",
    chains = 1,
    iter = 20,
    warmup = 10,
    refresh = 0,
    seed = 123L
  )
  expect_true(S7::S7_inherits(fit_sample, IdealAdaptorStanfit))
  expect_equal(fit_sample@backend, "cmdstanr")
  expect_true(inherits(fit_sample@stanfit, "stanfit"))
  expect_equal(nrow(rstan::extract(fit_sample@stanfit)$kappa_0), 10)
})
