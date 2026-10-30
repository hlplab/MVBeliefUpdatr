test_that("deprecated wrappers forward to the S7 fit-input constructor", {
  exposure <- data.frame(
    category = factor(c("A", "B")),
    group = factor(c("g1", "g1")),
    cue1 = c(0.1, 0.4),
    cue2 = c(0.2, 0.3)
  )
  test <- data.frame(
    response_category = factor(c("A", "B")),
    group = factor(c("g1", "g1")),
    cue1 = c(0.15, 0.45),
    cue2 = c(0.25, 0.35)
  )

  expect_warning(
    obj <- make_staninput(
      exposure = exposure,
      test = test,
      cues = c("cue1", "cue2"),
      category = "category",
      response = "response_category",
      group = "group",
      stanmodel = "NIW_ideal_adaptor"
    ),
    "new_ideal_adaptor_stanfit_input"
  )

  expect_true(S7::S7_inherits(obj, IdealAdaptorStanfitInput))
  expect_true(
    S7::S7_inherits(.new_staninput(obj), IdealAdaptorStaninput)
  )
})

test_that(
  "legacy wrappers trigger deprecation warnings for fixed parameters",
  {
  exposure <- data.frame(
    category = factor(c("A", "B")),
    group = factor(c("g1", "g1")),
    cue1 = c(0, 1),
    cue2 = c(0, 1)
  )
  test <- data.frame(
    response_category = factor(c("A", "B")),
    group = factor(c("g1", "g1")),
    cue1 = c(0.1, 0.9),
    cue2 = c(0.1, 0.9)
  )

  # lapse_rate triggers deprecation warning
  expect_warning(
    make_ideal_adaptor_stanfit_input(
      exposure = exposure,
      test = test,
      cues = c("cue1", "cue2"),
      category = "category",
      response = "response_category",
      group = "group",
      lapse_rate = 0.05,
      control = control_staninput(transform_type = "identity"),
      stanmodel = "NIW_ideal_adaptor"
    ),
    "lapse_rate"
  )

  # make_staninput with mu_0 and Sigma_0 triggers deprecation warnings
  mu_0_list <- list(A = c(0, 0), B = c(1, 1))
  Sigma_0_list <- list(A = diag(2), B = diag(2))
  expect_warning(
    make_staninput(
      exposure = exposure,
      test = test,
      cues = c("cue1", "cue2"),
      category = "category",
      response = "response_category",
      group = "group",
      mu_0 = mu_0_list,
      stanmodel = "NIW_ideal_adaptor"
    ),
    "mu_0"
  )
  expect_warning(
    make_staninput(
      exposure = exposure,
      test = test,
      cues = c("cue1", "cue2"),
      category = "category",
      response = "response_category",
      group = "group",
      Sigma_0 = Sigma_0_list,
      stanmodel = "NIW_ideal_adaptor"
    ),
    "Sigma_0"
  )
})

test_that(
  "legacy wrappers trigger deprecation warning for tau_scale and L_omega_eta",
  {
  exposure <- data.frame(
    category = factor(c("A", "B")),
    group = factor(c("g1", "g1")),
    cue1 = c(0.1, 0.4),
    cue2 = c(0.2, 0.3)
  )
  test <- data.frame(
    response_category = factor(c("A", "B")),
    group = factor(c("g1", "g1")),
    cue1 = c(0.15, 0.45),
    cue2 = c(0.25, 0.35)
  )

  expect_warning(
    make_staninput(
      exposure = exposure,
      test = test,
      cues = c("cue1", "cue2"),
      category = "category",
      response = "response_category",
      group = "group",
      tau_scale = 5
    ),
    "tau_scale"
  )

  expect_warning(
    make_ideal_adaptor_stanfit_input(
      exposure = exposure,
      test = test,
      cues = c("cue1", "cue2"),
      category = "category",
      response = "response_category",
      group = "group",
      L_omega_eta = 2
    ),
    "L_omega_eta"
  )
})

