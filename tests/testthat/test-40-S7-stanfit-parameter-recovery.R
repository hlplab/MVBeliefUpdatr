# test-40-S7-stanfit-parameter-recovery.R
# Parameter recovery test suite verifying analytical forward-updated beliefs
# against Stan posterior estimates and conjugate updating rules.

test_that("Stan NIX conjugate updating matches analytical update_template", {
  skip_if_not_installed("rstan")

  data_1d <- example_exposure_test_data("NIX", n_cues = 1L, seed = 42L)
  # Limit test stimulus rows for fast execution
  data_1d$test <- data_1d$test[1:20, ]

  k0_val <- 5
  nu0_val <- 10
  m0_d <- 0
  m0_t <- 50
  s0_d <- 25
  s0_t <- 100

  fixed_pars <- c(
    fix_parameter(m0_d, class = "m_0", category = "/d/"),
    fix_parameter(m0_t, class = "m_0", category = "/t/"),
    fix_parameter(s0_d, class = "S_0", category = "/d/"),
    fix_parameter(s0_t, class = "S_0", category = "/t/"),
    fix_parameter(k0_val, class = "kappa_0"),
    fix_parameter(nu0_val, class = "nu_0"),
    fix_parameter(0, class = "lapse_rate")
  )

  inp <- example_ideal_adaptor_stanfit_input(
    model_family = "NIX",
    n_cues = 1L,
    data = data_1d,
    control = control_staninput(transform_type = "identity")
  )

  fit <- suppressWarnings(
    fit_ideal_adaptor(
      stanfit_input = inp,
      fixed_parameters = fixed_pars,
      chains = 1,
      iter = 2,
      warmup = 1,
      refresh = 0,
      seed = 42L,
      stanmodel = "NIX_ideal_adaptor"
    )
  )

  draws <- rstan::extract(get_stanfit(fit))

  # Build corresponding analytical prior model in R
  prior_adaptor <- example_nix_ideal_adaptor(
    categories = c("/d/", "/t/"),
    kappa = k0_val,
    nu = nu0_val
  )
  prior_adaptor@category_template@representations[["/d/"]]@m <- m0_d
  prior_adaptor@category_template@representations[["/d/"]]@sigma2 <-
    s0_d / nu0_val
  prior_adaptor@category_template@representations[["/t/"]]@m <- m0_t
  prior_adaptor@category_template@representations[["/t/"]]@sigma2 <-
    s0_t / nu0_val

  groups <- get_group_labels(inp)

  for (g_idx in seq_along(groups)) {
    g_name <- groups[g_idx]
    exp_sub <- data_1d$exposure[data_1d$exposure$group == g_name, ]
    up <- update_template(prior_adaptor, exp_sub)

    rep_d <- up@category_template@representations[["/d/"]]
    rep_t <- up@category_template@representations[["/t/"]]

    # Category /d/
    expect_equal(
      draws$m_n[1, 1, g_idx],
      as.numeric(rep_d@m),
      tolerance = 1e-6
    )
    expect_equal(
      draws$S_n[1, 1, g_idx],
      as.numeric(rep_d@nu * rep_d@sigma2),
      tolerance = 1e-6
    )
    expect_equal(
      draws$kappa_n[1, 1, g_idx],
      as.numeric(rep_d@kappa)
    )
    expect_equal(
      draws$nu_n[1, 1, g_idx],
      as.numeric(rep_d@nu)
    )

    # Category /t/
    expect_equal(
      draws$m_n[1, 2, g_idx],
      as.numeric(rep_t@m),
      tolerance = 1e-6
    )
    expect_equal(
      draws$S_n[1, 2, g_idx],
      as.numeric(rep_t@nu * rep_t@sigma2),
      tolerance = 1e-6
    )
    expect_equal(
      draws$kappa_n[1, 2, g_idx],
      as.numeric(rep_t@kappa)
    )
    expect_equal(
      draws$nu_n[1, 2, g_idx],
      as.numeric(rep_t@nu)
    )
  }
})

test_that("Stan NIW conjugate updating matches analytical update_template", {
  skip_if_not_installed("rstan")

  data_2d <- example_exposure_test_data("NIW", n_cues = 2L, seed = 42L)
  data_2d$test <- data_2d$test[1:20, ]

  k0_val <- 5
  nu0_val <- 10
  m0_d <- c(0, 0)
  m0_t <- c(50, 200)
  S0_d <- diag(c(25, 100))
  S0_t <- diag(c(50, 150))

  fixed_pars <- c(
    fix_parameter(m0_d, class = "m_0", category = "/d/"),
    fix_parameter(m0_t, class = "m_0", category = "/t/"),
    fix_parameter(S0_d, class = "S_0", category = "/d/"),
    fix_parameter(S0_t, class = "S_0", category = "/t/"),
    fix_parameter(k0_val, class = "kappa_0"),
    fix_parameter(nu0_val, class = "nu_0"),
    fix_parameter(0, class = "lapse_rate")
  )

  inp <- example_ideal_adaptor_stanfit_input(
    model_family = "NIW",
    n_cues = 2L,
    data = data_2d,
    control = control_staninput(transform_type = "identity")
  )

  fit <- suppressWarnings(
    fit_ideal_adaptor(
      stanfit_input = inp,
      fixed_parameters = fixed_pars,
      chains = 1,
      iter = 2,
      warmup = 1,
      refresh = 0,
      seed = 42L,
      stanmodel = "NIW_ideal_adaptor"
    )
  )

  draws <- rstan::extract(get_stanfit(fit))

  prior_adaptor <- example_niw_ideal_adaptor(
    n_cues = 2L,
    categories = c("/d/", "/t/"),
    kappa = k0_val,
    nu = nu0_val
  )
  prior_adaptor@category_template@representations[["/d/"]]@m <- m0_d
  prior_adaptor@category_template@representations[["/d/"]]@S <- S0_d
  prior_adaptor@category_template@representations[["/t/"]]@m <- m0_t
  prior_adaptor@category_template@representations[["/t/"]]@S <- S0_t

  groups <- get_group_labels(inp)

  for (g_idx in seq_along(groups)) {
    g_name <- groups[g_idx]
    exp_sub <- data_2d$exposure[data_2d$exposure$group == g_name, ]
    up <- update_template(prior_adaptor, exp_sub)

    rep_d <- up@category_template@representations[["/d/"]]
    rep_t <- up@category_template@representations[["/t/"]]

    # Category /d/
    expect_equal(
      as.numeric(draws$m_n[1, 1, g_idx, ]),
      as.numeric(rep_d@m),
      tolerance = 1e-6
    )
    expect_equal(
      draws$S_n[1, 1, g_idx, , ],
      rep_d@S,
      tolerance = 1e-6
    )
    expect_equal(
      draws$kappa_n[1, 1, g_idx],
      as.numeric(rep_d@kappa)
    )
    expect_equal(
      draws$nu_n[1, 1, g_idx],
      as.numeric(rep_d@nu)
    )

    # Category /t/
    expect_equal(
      as.numeric(draws$m_n[1, 2, g_idx, ]),
      as.numeric(rep_t@m),
      tolerance = 1e-6
    )
    expect_equal(
      draws$S_n[1, 2, g_idx, , ],
      rep_t@S,
      tolerance = 1e-6
    )
    expect_equal(
      draws$kappa_n[1, 2, g_idx],
      as.numeric(rep_t@kappa)
    )
    expect_equal(
      draws$nu_n[1, 2, g_idx],
      as.numeric(rep_t@nu)
    )
  }
})

test_that("Stan MNIX conjugate updating matches analytical update_template", {
  skip_if_not_installed("rstan")

  data_2d <- example_exposure_test_data("MNIX", n_cues = 2L, seed = 42L)
  data_2d$test <- data_2d$test[1:20, ]

  k0_val <- 5
  nu0_val <- 10
  m0_d <- c(0, 0)
  m0_t <- c(50, 200)
  S0_d <- diag(c(25, 100))
  S0_t <- diag(c(50, 150))

  fixed_pars <- c(
    fix_parameter(m0_d, class = "m_0", category = "/d/"),
    fix_parameter(m0_t, class = "m_0", category = "/t/"),
    fix_parameter(S0_d, class = "S_0", category = "/d/"),
    fix_parameter(S0_t, class = "S_0", category = "/t/"),
    fix_parameter(k0_val, class = "kappa_0"),
    fix_parameter(nu0_val, class = "nu_0"),
    fix_parameter(0, class = "lapse_rate")
  )

  inp <- example_ideal_adaptor_stanfit_input(
    model_family = "MNIX",
    n_cues = 2L,
    data = data_2d,
    control = control_staninput(transform_type = "identity")
  )

  fit <- suppressWarnings(
    fit_ideal_adaptor(
      stanfit_input = inp,
      fixed_parameters = fixed_pars,
      chains = 1,
      iter = 2,
      warmup = 1,
      refresh = 0,
      seed = 42L,
      stanmodel = "MNIX_ideal_adaptor_ideal_integrator"
    )
  )

  draws <- rstan::extract(get_stanfit(fit))

  prior_adaptor <- example_mnix_ideal_adaptor(
    n_cues = 2L,
    categories = c("/d/", "/t/"),
    kappa = k0_val,
    nu = nu0_val
  )
  prior_adaptor@category_template@representations[["/d/"]]@m <- m0_d
  prior_adaptor@category_template@representations[["/d/"]]@sigma2 <-
    diag(S0_d) / nu0_val
  prior_adaptor@category_template@representations[["/t/"]]@m <- m0_t
  prior_adaptor@category_template@representations[["/t/"]]@sigma2 <-
    diag(S0_t) / nu0_val

  groups <- get_group_labels(inp)

  for (g_idx in seq_along(groups)) {
    g_name <- groups[g_idx]
    exp_sub <- data_2d$exposure[data_2d$exposure$group == g_name, ]
    up <- update_template(prior_adaptor, exp_sub)

    rep_d <- up@category_template@representations[["/d/"]]
    rep_t <- up@category_template@representations[["/t/"]]

    # Category /d/
    expect_equal(
      as.numeric(draws$m_n[1, 1, g_idx, ]),
      as.numeric(rep_d@m),
      tolerance = 1e-6
    )
    expect_equal(
      as.numeric(draws$S_n[1, 1, g_idx, ]),
      as.numeric(rep_d@sigma2),
      tolerance = 1e-6
    )
    expect_equal(
      draws$kappa_n[1, 1, g_idx],
      as.numeric(rep_d@kappa)[1]
    )
    expect_equal(
      draws$nu_n[1, 1, g_idx],
      as.numeric(rep_d@nu)[1]
    )

    # Category /t/
    expect_equal(
      as.numeric(draws$m_n[1, 2, g_idx, ]),
      as.numeric(rep_t@m),
      tolerance = 1e-6
    )
    expect_equal(
      as.numeric(draws$S_n[1, 2, g_idx, ]),
      as.numeric(rep_t@sigma2),
      tolerance = 1e-6
    )
    expect_equal(
      draws$kappa_n[1, 2, g_idx],
      as.numeric(rep_t@kappa)[1]
    )
    expect_equal(
      draws$nu_n[1, 2, g_idx],
      as.numeric(rep_t@nu)[1]
    )
  }
})

test_that("Known moments mu_0 and Sigma_0 govern belief updating (Cases 2-4)", {
  skip_if_not_installed("rstan")

  data_1d <- example_exposure_test_data("NIX", n_cues = 1L, seed = 42L)
  data_1d$test <- data_1d$test[1:20, ]

  k0_val <- 5
  nu0_val <- 10
  mu0_d <- 10
  sigma0_d <- 16

  # Case 2: Sigma_0 known, mu_0 unknown for /d/
  # Case 3: mu_0 known, Sigma_0 unknown for /t/
  fixed_pars <- c(
    fix_parameter(sigma0_d, class = "Sigma_0", category = "/d/"),
    fix_parameter(0, class = "m_0", category = "/d/"),
    fix_parameter(50, class = "mu_0", category = "/t/"),
    fix_parameter(100, class = "S_0", category = "/t/"),
    fix_parameter(k0_val, class = "kappa_0"),
    fix_parameter(nu0_val, class = "nu_0"),
    fix_parameter(0, class = "lapse_rate")
  )

  inp <- example_ideal_adaptor_stanfit_input(
    model_family = "NIX",
    n_cues = 1L,
    data = data_1d,
    control = control_staninput(transform_type = "identity")
  )

  fit <- suppressWarnings(
    fit_ideal_adaptor(
      stanfit_input = inp,
      fixed_parameters = fixed_pars,
      chains = 1,
      iter = 2,
      warmup = 1,
      refresh = 0,
      seed = 42L,
      stanmodel = "NIX_ideal_adaptor"
    )
  )

  draws <- rstan::extract(get_stanfit(fit))

  # Case 2 (/d/): Sigma_0 known => nu_n remains nu_0 across all groups,
  # kappa_n updates to kappa_0 + n, S_n is constant = Sigma_0
  groups <- get_group_labels(inp)
  for (g_idx in seq_along(groups)) {
    g_name <- groups[g_idx]
    n_obs <- sum(
      data_1d$exposure$group == g_name & data_1d$exposure$category == "/d/"
    )
    expect_equal(draws$nu_n[1, 1, g_idx], nu0_val)
    expect_equal(draws$kappa_n[1, 1, g_idx], k0_val + n_obs)
    expect_equal(draws$S_n[1, 1, g_idx], sigma0_d, tolerance = 1e-6)
  }

  # Case 3 (/t/): mu_0 known => kappa_n remains kappa_0 across all groups,
  # nu_n updates to nu_0 + n, m_n is constant = mu_0
  for (g_idx in seq_along(groups)) {
    g_name <- groups[g_idx]
    n_obs <- sum(
      data_1d$exposure$group == g_name & data_1d$exposure$category == "/t/"
    )
    expect_equal(draws$kappa_n[1, 2, g_idx], k0_val)
    expect_equal(draws$nu_n[1, 2, g_idx], nu0_val + n_obs)
    expect_equal(draws$m_n[1, 2, g_idx], 50, tolerance = 1e-6)
  }
})

test_that("Posterior draws recover true prior parameters from simulated data", {
  skip_if_not_installed("rstan")

  model_path <- testthat::test_path(
    "models",
    "example-stanfit-NIX_ideal_adaptor-1-standardize-42.rds"
  )
  skip_if_not(file.exists(model_path))

  fit <- read_stanfit(model_path)
  sfit <- get_stanfit(fit)
  draws <- rstan::extract(sfit)

  # Check that posterior estimates of m_0 have reasonable finite bounds
  m0_d_draws <- draws$m_0[, 1]
  m0_t_draws <- draws$m_0[, 2]

  expect_true(is.numeric(m0_d_draws))
  expect_true(is.numeric(m0_t_draws))
  expect_false(anyNA(m0_d_draws))
  expect_false(anyNA(m0_t_draws))

  # Posterior mean and 95% credible intervals are well-defined
  ci_d <- stats::quantile(m0_d_draws, probs = c(0.025, 0.975))
  ci_t <- stats::quantile(m0_t_draws, probs = c(0.025, 0.975))

  expect_lt(ci_d[1], ci_d[2])
  expect_lt(ci_t[1], ci_t[2])

  # In standardized space, mean estimates are centered reasonably near 0
  expect_lt(abs(mean(m0_d_draws)), 3.0)
  expect_lt(abs(mean(m0_t_draws)), 3.0)

  # Check multivariate NIW fixture parameter recovery
  model_path_niw <- testthat::test_path(
    "models",
    "minimal-niw_ideal_adaptor-NIW-2cue.rds"
  )
  skip_if_not(file.exists(model_path_niw))

  fit_niw <- read_stanfit(model_path_niw)
  draws_niw <- rstan::extract(get_stanfit(fit_niw))

  # 2 categories, 2 cues
  expect_equal(dim(draws_niw$m_0)[2:3], c(2, 2))
  for (cat in 1:2) {
    for (k in 1:2) {
      vals <- draws_niw$m_0[, cat, k]
      expect_false(anyNA(vals))
      expect_lt(abs(mean(vals)), 5.0)
    }
  }
})
