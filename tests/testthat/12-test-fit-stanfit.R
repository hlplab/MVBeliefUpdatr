make_minimal_stanfit_data <- function(cues = c("VOT")) {
  model <- example_MVG_ideal_observer(
    example = 1,
    categories = c("/b/", "/p/"),
    cues = cues
  )
  data <- sample_data_from_model(model, Ns = c(10, 10), randomize.order = FALSE)

  exposure <- data[1:10, , drop = FALSE]
  test <- data[11:20, , drop = FALSE]

  exposure$group <- factor(rep(c("g1", "g2"), each = 5), levels = c("g1", "g2"))
  exposure$category <- factor(exposure$category, levels = c("/b/", "/p/"))
  exposure$response <- factor(exposure$category, levels = c("/b/", "/p/"))

  test$group <- factor(rep(c("g1", "g2"), each = 5), levels = c("g1", "g2"))
  test$category <- factor(test$category, levels = c("/b/", "/p/"))
  test$response <- factor(test$category, levels = c("/b/", "/p/"))

  list(exposure = exposure, test = test)
}

test_that("legacy fit path produces an ideal-adaptor stanfit input for minimal data", {
  skip_if_not_installed("rstan")

  data <- make_minimal_stanfit_data(cues = "VOT")
  staninput <- new_ideal_adaptor_staninput(
    exposure = data$exposure,
    test = data$test,
    cues = "VOT",
    category = "category",
    response = "response",
    group = "group",
    control = control_staninput(transform_type = "identity"),
    stanmodel = "NIW_ideal_adaptor"
  )

  expect_true(is.ideal_adaptor_stanfit_input(staninput))
  expect_true(is.data.frame(staninput$data))
  expect_true(is.list(staninput$staninput$transformed))
  expect_true(is.list(staninput$staninput$untransformed))
})

test_that("legacy fit path supports the new constructor API for NIX and MNIX minimal data", {
  skip_if_not_installed("rstan")

  nix_data <- make_minimal_stanfit_data(cues = "VOT")
  nix_staninput <- new_ideal_adaptor_staninput(
    exposure = nix_data$exposure,
    test = nix_data$test,
    cues = "VOT",
    category = "category",
    response = "response",
    group = "group",
    control = control_staninput(transform_type = "identity"),
    stanmodel = "NIX_ideal_adaptor"
  )

  expect_true(is.ideal_adaptor_stanfit_input(nix_staninput))
  expect_true(all(c("x_mean_exposure", "x_sd_exposure") %in% names(nix_staninput$staninput$untransformed)))

  mnix_data <- make_minimal_stanfit_data(cues = c("VOT", "f0_semitones"))
  mnix_staninput <- new_ideal_adaptor_staninput(
    exposure = mnix_data$exposure,
    test = mnix_data$test,
    cues = c("VOT", "f0_semitones"),
    category = "category",
    response = "response",
    group = "group",
    control = control_staninput(transform_type = "identity"),
    stanmodel = "MNIX_ideal_adaptor"
  )

  expect_true(is.ideal_adaptor_stanfit_input(mnix_staninput))
  expect_true(all(c("x_mean_exposure", "x_ss_exposure") %in% names(mnix_staninput$staninput$untransformed)))
})
