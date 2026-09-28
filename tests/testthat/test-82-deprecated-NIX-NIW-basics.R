test_that("get_NIW_posterior_predictive calculates correct densities and warns deprecation", {
  m <- c(0, 0)
  S <- diag(2)
  kappa <- 5
  nu <- 10

  x <- matrix(c(0, 0, 1, 1), nrow = 2, byrow = TRUE)
  expect_warning(
    log_dens <- get_NIW_posterior_predictive(x, m, S, kappa, nu, log = TRUE),
    class = "lifecycle_warning_deprecated"
  )
  expect_length(log_dens, 2)
  expect_true(all(is.finite(log_dens)))
  expect_true(log_dens[1] > log_dens[2]) # (0,0) is closer to mean

  expect_warning(
    dens <- get_NIW_posterior_predictive(x, m, S, kappa, nu, log = FALSE),
    class = "lifecycle_warning_deprecated"
  )
  expect_equal(dens, exp(log_dens))

  # Test noise treatment
  Sigma_noise <- diag(c(0.5, 0.5))
  expect_warning(
    log_dens_noise <- get_NIW_posterior_predictive(
      x, m, S, kappa, nu,
      Sigma_noise = Sigma_noise,
      noise_treatment = "marginalize",
      log = TRUE
    ),
    class = "lifecycle_warning_deprecated"
  )
  expect_length(log_dens_noise, 2)
  expect_true(all(is.finite(log_dens_noise)))
})

test_that("get_NIX_posterior_predictive calculates correct 1D densities and warns deprecation", {
  m <- 0
  sigma2 <- 1
  kappa <- 5
  nu <- 10

  x <- c(0, 1, 2)
  expect_warning(
    log_dens <- get_NIX_posterior_predictive(x, m, sigma2, kappa, nu, log = TRUE),
    class = "lifecycle_warning_deprecated"
  )
  expect_length(log_dens, 3)
  expect_true(all(is.finite(log_dens)))
  expect_true(log_dens[1] > log_dens[2])
  expect_true(log_dens[2] > log_dens[3])

  expect_warning(
    dens <- get_NIX_posterior_predictive(x, m, sigma2, kappa, nu, log = FALSE),
    class = "lifecycle_warning_deprecated"
  )
  expect_equal(dens, exp(log_dens))

  # Test noise treatment
  expect_warning(
    log_dens_noise <- get_NIX_posterior_predictive(
      x, m, sigma2, kappa, nu,
      Sigma_noise = 0.5,
      noise_treatment = "marginalize",
      log = TRUE
    ),
    class = "lifecycle_warning_deprecated"
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
