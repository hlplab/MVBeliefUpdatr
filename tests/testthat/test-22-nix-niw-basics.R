test_that("get_expected_mu_from_m and get_m_from_expected_mu work correctly", {
  m_val <- c(0, 1)
  expect_equal(get_expected_mu_from_m(m_val), m_val)
  expect_equal(get_m_from_expected_mu(m_val), m_val)

  m_list <- list(c(0, 1), c(2, 3))
  expect_equal(get_expected_mu_from_m(m_list), m_list)
  expect_equal(get_m_from_expected_mu(m_list), m_list)
})

test_that("get_expected_Sigma_from_S and get_S_from_expected_Sigma work correctly", {
  S_mat <- diag(c(4, 9))
  nu_val <- 10
  D <- 2
  # Sigma = S / (nu - D - 1) = S / 7
  expected_Sigma <- S_mat / 7
  expect_equal(get_expected_Sigma_from_S(S_mat, nu_val), expected_Sigma)
  expect_equal(get_S_from_expected_Sigma(expected_Sigma, nu_val), S_mat)
})

test_that("get_NIW_posterior_predictive calculates correct densities", {
  m <- c(0, 0)
  S <- diag(2)
  kappa <- 5
  nu <- 10

  x <- matrix(c(0, 0, 1, 1), nrow = 2, byrow = TRUE)
  log_dens <- get_NIW_posterior_predictive(x, m, S, kappa, nu, log = TRUE)
  expect_length(log_dens, 2)
  expect_true(all(is.finite(log_dens)))
  expect_true(log_dens[1] > log_dens[2]) # (0,0) is closer to mean

  dens <- get_NIW_posterior_predictive(x, m, S, kappa, nu, log = FALSE)
  expect_equal(dens, exp(log_dens))

  # Test noise treatment
  Sigma_noise <- diag(c(0.5, 0.5))
  log_dens_noise <- get_NIW_posterior_predictive(
    x, m, S, kappa, nu,
    Sigma_noise = Sigma_noise,
    noise_treatment = "marginalize",
    log = TRUE
  )
  expect_length(log_dens_noise, 2)
  expect_true(all(is.finite(log_dens_noise)))
})

test_that("get_NIX_posterior_predictive calculates correct 1D densities", {
  m <- 0
  sigma2 <- 1
  kappa <- 5
  nu <- 10

  x <- c(0, 1, 2)
  log_dens <- get_NIX_posterior_predictive(x, m, sigma2, kappa, nu, log = TRUE)
  expect_length(log_dens, 3)
  expect_true(all(is.finite(log_dens)))
  expect_true(log_dens[1] > log_dens[2])
  expect_true(log_dens[2] > log_dens[3])

  dens <- get_NIX_posterior_predictive(x, m, sigma2, kappa, nu, log = FALSE)
  expect_equal(dens, exp(log_dens))

  # Test noise treatment
  log_dens_noise <- get_NIX_posterior_predictive(
    x, m, sigma2, kappa, nu,
    Sigma_noise = 0.5,
    noise_treatment = "marginalize",
    log = TRUE
  )
  expect_length(log_dens_noise, 3)
  expect_true(all(is.finite(log_dens_noise)))
})

test_that("get_D emits deprecation warning and returns length of cue labels", {
  model <- example_niw_ideal_adaptor(n_cues = 2)
  expect_warning(
    d_val <- get_D(model),
    class = "lifecycle_warning_deprecated"
  )
  expect_equal(d_val, 2)
})
