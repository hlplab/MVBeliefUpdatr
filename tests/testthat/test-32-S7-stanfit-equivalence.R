test_that(
  "1-cue NIX and NIW stanfits yield equivalent inputs and posteriors",
  {
    skip_if_not_installed("rstan")

    # Generate identical 1D dataset
    data_1d <- example_exposure_test_data(
      model_family = "NIX",
      n_cues = 1L,
      seed = 42L
    )

    input_nix <- example_ideal_adaptor_stanfit_input(
      model_family = "NIX",
      n_cues = 1L,
      data = data_1d,
      control = control_staninput(transform_type = "identity")
    )

    input_niw <- example_ideal_adaptor_stanfit_input(
      model_family = "NIW",
      n_cues = 1L,
      data = data_1d,
      control = control_staninput(transform_type = "identity")
    )

  # Check staninput equivalence for exposure sufficient stats
  expect_equal(get_cue_labels(input_nix), get_cue_labels(input_niw))
  expect_equal(get_category_labels(input_nix), get_category_labels(input_niw))
  expect_equal(get_group_labels(input_nix), get_group_labels(input_niw))
  expect_equal(get_labels(input_nix), get_labels(input_niw))

  # Exposure centered sum of squares equivalence
  expect_equal(
    as.vector(input_nix@staninput@values$x_ss_exposure[, 1]),
    as.vector(input_niw@staninput@values$x_ss_exposure[, 1, 1, 1]),
    tolerance = 1e-6
  )

  # Check exact mathematical log posterior density equivalence via
  # rstan::log_prob
  mod_nix <- stanmodels$NIX_ideal_adaptor
  mod_niw <- stanmodels$NIW_ideal_adaptor

  fit_nix_init <- rstan::sampling(
    mod_nix,
    data = input_nix@staninput@values,
    iter = 1,
    chains = 1,
    warmup = 0,
    refresh = 0
  )
  fit_niw_init <- rstan::sampling(
    mod_niw,
    data = input_niw@staninput@values,
    iter = 1,
    chains = 1,
    warmup = 0,
    refresh = 0
  )

  init_params_nix <- list(
    kappa_0 = 10,
    nu_0 = 15,
    lapse_rate_param = array(0.05, dim = 1),
    m_0_param = c(-0.5, 0.5),
    m_0_tau = array(1.5, dim = 1),
    tau_0_param = c(1.2, 1.8)
  )

  init_params_niw <- list(
    kappa_0 = 10,
    nu_0 = 15,
    lapse_rate_param = array(0.05, dim = 1),
    m_0_param = array(c(-0.5, 0.5), dim = c(2, 1)),
    m_0_tau = array(1.5, dim = 1),
    m_0_L_omega = matrix(1, nrow = 1, ncol = 1),
    tau_0_param = array(c(1.2, 1.8), dim = c(2, 1)),
    L_omega_0_param = array(1, dim = c(2, 1, 1))
  )

  u_nix <- rstan::unconstrain_pars(fit_nix_init, init_params_nix)
  u_niw <- rstan::unconstrain_pars(fit_niw_init, init_params_niw)

  expect_equal(u_nix, u_niw, tolerance = 1e-7)

  lp_nix <- rstan::log_prob(fit_nix_init, u_nix, adjust_transform = FALSE)
  lp_niw <- rstan::log_prob(fit_niw_init, u_niw, adjust_transform = FALSE)
  expect_equal(lp_nix, lp_niw, tolerance = 1e-6)

  lp_nix_adj <- rstan::log_prob(fit_nix_init, u_nix, adjust_transform = TRUE)
  lp_niw_adj <- rstan::log_prob(fit_niw_init, u_niw, adjust_transform = TRUE)
  expect_equal(lp_nix_adj, lp_niw_adj, tolerance = 1e-6)
})
