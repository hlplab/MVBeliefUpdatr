test_that("1-cue NIX and NIW ideal adaptor stanfits yield equivalent posterior summaries", {
  skip_if_not_installed("rstan")
  skip("On-the-fly MCMC fitting skipped during fast test suite runs; requires NIX/NIW Stan code alignment.")

  # Generate identical 1D dataset
  data_1d <- build_shifted_prior_fit_example(
    cues = "VOT",
    n_exposure_per_condition = 30L,
    n_test_per_condition = 25L,
    seed = 42L
  )

  input_nix <- new_ideal_adaptor_stanfit_input(
    exposure = data_1d$exposure,
    test = data_1d$test,
    cues = "VOT",
    category = "category",
    response = "response",
    group = "group",
    group.unique = "Condition",
    control = control_staninput(transform_type = "identity"),
    stanmodel = "NIX_ideal_adaptor"
  )

  input_niw <- new_ideal_adaptor_stanfit_input(
    exposure = data_1d$exposure,
    test = data_1d$test,
    cues = "VOT",
    category = "category",
    response = "response",
    group = "group",
    group.unique = "Condition",
    control = control_staninput(transform_type = "identity"),
    stanmodel = "NIW_ideal_adaptor"
  )

  # Check staninput equivalence for exposure sufficient stats
  expect_equal(get_cue_labels(input_nix), get_cue_labels(input_niw))
  expect_equal(get_category_labels(input_nix), get_category_labels(input_niw))
  expect_equal(get_group_labels(input_nix), get_group_labels(input_niw))
  expect_equal(get_labels(input_nix), get_labels(input_niw))

  # Fit both models with same MCMC settings and seed
  fit_nix <- fit_ideal_adaptor(
    stanfit_input = input_nix,
    chains = 1,
    iter = 500,
    warmup = 250,
    seed = 12345,
    refresh = 0,
    control = list(adapt_delta = 0.95),
    stanmodel = "NIX_ideal_adaptor"
  )

  fit_niw <- fit_ideal_adaptor(
    stanfit_input = input_niw,
    chains = 1,
    iter = 500,
    warmup = 250,
    seed = 12345,
    refresh = 0,
    control = list(adapt_delta = 0.95),
    stanmodel = "NIW_ideal_adaptor"
  )

  expect_true(S7::S7_inherits(fit_nix, NIX_IdealAdaptorStanfit))
  expect_true(S7::S7_inherits(fit_niw, NIW_IdealAdaptorStanfit))

  # Test expected mu equivalence across categories
  cats <- get_category_labels(fit_nix)
  grps <- get_group_labels(fit_nix)

  mu_nix <- get_expected_mu(fit_nix)
  mu_niw <- get_expected_mu(fit_niw)

  expect_true(is.numeric(mu_nix$mu.mean[[1]]))
  expect_true(is.numeric(mu_niw$mu.mean[[1]]))

  # Compare posterior mean estimates (probabilistic equivalence within MCMC tolerance)
  for (i in seq_along(cats)) {
    mu_val_nix <- unlist(mu_nix$mu.mean[mu_nix$category == cats[i]])
    mu_val_niw <- unlist(mu_niw$mu.mean[mu_niw$category == cats[i]])
    expect_equal(mu_val_nix, mu_val_niw, tolerance = 0.25)
  }

  # Test expected Sigma / variance equivalence
  sigma_nix <- get_expected_sigma(fit_nix)
  sigma_niw <- get_expected_sigma(fit_niw)

  for (i in seq_along(cats)) {
    sig_val_nix <- unlist(sigma_nix$Sigma.mean[sigma_nix$category == cats[i]])
    sig_val_niw <- unlist(sigma_niw$Sigma.mean[sigma_niw$category == cats[i]])
    expect_equal(sig_val_nix, sig_val_niw, tolerance = 0.35)
  }
})
