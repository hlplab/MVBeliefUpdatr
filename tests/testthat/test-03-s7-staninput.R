# Build a small, general-purpose synthetic exposure/test dataset for broad
# constructor and Stan-input regression tests. The helper can be configured to
# vary the number of exposure observations, the number of test observations,
# the number of groups, and the number of categories.
make_minimal_staninput_data <- function(
  cues = c("cue1"),
  n_obs_exposure = 8L,
  n_obs_test = 8L,
  n_group = 2L,
  n_category = 2L
) {
  cue_names <- as.character(cues)
  category_levels <- paste0("A", seq_len(n_category))
  group_levels <- paste0("g", seq_len(n_group))

  exposure <- data.frame(
    group = factor(rep(group_levels, length.out = n_obs_exposure), levels = group_levels),
    category = factor(rep(category_levels, length.out = n_obs_exposure), levels = category_levels),
    response = factor(rep(category_levels, length.out = n_obs_exposure), levels = category_levels)
  )
  if (n_obs_exposure == 0L) {
    exposure <- exposure[0, , drop = FALSE]
  }

  if (length(cue_names) > 0) {
    base_values <- c(0.12, 0.82, 0.24, 0.76, 0.36, 0.64, 0.48, 0.88)
    for (i in seq_along(cue_names)) {
      cue_name <- cue_names[[i]]
      if (n_obs_exposure > 0L) {
        cue_values <- pmin(0.98, rep(base_values, length.out = n_obs_exposure) + 0.03 * (i - 1))
      } else {
        cue_values <- numeric(0)
      }
      exposure[[cue_name]] <- cue_values
    }
  }

  test <- data.frame(
    group = factor(rep(group_levels, length.out = n_obs_test), levels = group_levels),
    response = factor(rep(category_levels, length.out = n_obs_test), levels = category_levels)
  )
  if (n_obs_test == 0L) {
    test <- test[0, , drop = FALSE]
  }

  test_values <- c(0.15, 0.85, 0.30, 0.70, 0.45, 0.55, 0.60, 0.90)
  for (i in seq_along(cue_names)) {
    cue_name <- cue_names[[i]]
    if (n_obs_test > 0L) {
      test[[cue_name]] <- pmin(0.98, rep(test_values, length.out = n_obs_test) + 0.01 * (i - 1))
    } else {
      test[[cue_name]] <- numeric(0)
    }
  }

  list(exposure = exposure, test = test)
}

test_that("new_ideal_adaptor_staninput accepts empty exposure data", {
  exposure <- data.frame(
    category = factor(character(), levels = "A"),
    group = factor(character(), levels = "g1"),
    cue1 = numeric(),
    cue2 = numeric()
  )
  test <- data.frame(
    response = factor("A", levels = "A"),
    group = factor("g1", levels = "g1"),
    cue1 = 1,
    cue2 = 2
  )

  control <- list(
    tau_scale = 1,
    L_omega_eta = 1,
    split_loglik_per_observation = 0,
    transform_type = "identity"
  )

  res <- new_ideal_adaptor_staninput(
    exposure = exposure,
    test = test,
    cues = c("cue1", "cue2"),
    category = "category",
    response = "response",
    group = "group",
    control = control,
    stanmodel = "NIW_ideal_adaptor"
  )

  expect_true(is.list(res))
  expect_true(all(c("staninput", "data", "transform_information") %in% names(res)))
  expect_true(is.data.frame(res$data))
  expect_true(is.list(res$staninput))
  expect_true(all(c("transformed", "untransformed") %in% names(res$staninput)))
  expect_true(is.list(res$transform_information))
  expect_equal(res$staninput$transformed$N_test, 1)
})

test_that("new_ideal_adaptor_staninput preserves group.unique in exposure data when it's specified", {
  exposure <- data.frame(
    Condition = factor("baseline"),
    group = factor("g1"),
    category = factor("A", levels = "A"),
    cue1 = 1
  )
  test <- data.frame(
    Condition = factor("baseline"),
    group = factor("g1"),
    response = factor("A", levels = "A"),
    cue1 = 1.1
  )

  control <- list(
    tau_scale = 1,
    L_omega_eta = 1,
    split_loglik_per_observation = 0,
    transform_type = "identity"
  )

  res <- expect_no_error(
    new_ideal_adaptor_staninput(
      exposure = exposure,
      test = test,
      cues = "cue1",
      category = "category",
      response = "response",
      group = "group",
      group.unique = "Condition",
      control = control,
      stanmodel = "NIW_ideal_adaptor"
    )
  )

  expect_true("Condition" %in% names(res$data))
  expect_true(all(res$data$Condition[res$data$Phase == "exposure"] == "baseline"))
})

test_that("synthetic exposure data has non-zero variance", {
  data <- make_minimal_staninput_data(cues = c("cue1", "cue2"))

  expect_gt(stats::var(data$exposure$cue1), 0)
  expect_gt(stats::var(data$exposure$cue2), 0)
})

test_that("new_ideal_adaptor_staninput emits expected prior-array shapes for NIX, NIW, and MNIX", {
  data <- make_minimal_staninput_data(cues = c("cue1", "cue2"))

  nix_res <- new_ideal_adaptor_staninput(
    exposure = data$exposure,
    test = data$test,
    cues = "cue1",
    category = "category",
    response = "response",
    group = "group",
    control = control_staninput(transform_type = "identity"),
    stanmodel = "NIX_ideal_adaptor",
    fix_parameters = list(
      mu_0 = lapply(seq_len(2), function(i) 0.1),
      Sigma_0 = lapply(seq_len(2), function(i) 0.1)
    )
  )

  expect_equal(length(nix_res$staninput$transformed$mu_0_data), 2)
  expect_equal(length(nix_res$staninput$transformed$Sigma_0_data), 2)

  niw_res <- new_ideal_adaptor_staninput(
    exposure = data$exposure,
    test = data$test,
    cues = c("cue1", "cue2"),
    category = "category",
    response = "response",
    group = "group",
    control = control_staninput(transform_type = "identity"),
    stanmodel = "NIW_ideal_adaptor",
    fix_parameters = list(
      mu_0 = lapply(seq_len(2), function(i) c(0.1, 0.2)),
      Sigma_0 = lapply(seq_len(2), function(i) diag(2))
    )
  )

  expect_equal(dim(niw_res$staninput$transformed$mu_0_data), c(2, 2))
  expect_equal(dim(niw_res$staninput$transformed$Sigma_0_data), c(2, 2, 2))

  mnix_res <- new_ideal_adaptor_staninput(
    exposure = data$exposure,
    test = data$test,
    cues = c("cue1", "cue2"),
    category = "category",
    response = "response",
    group = "group",
    control = control_staninput(transform_type = "identity"),
    stanmodel = "MNIX_ideal_adaptor",
    fix_parameters = list(
      mu_0 = lapply(seq_len(2), function(i) c(0.1, 0.2)),
      Sigma_0 = lapply(seq_len(2), function(i) diag(2))
    )
  )

  expect_equal(dim(mnix_res$staninput$transformed$mu_0_data), c(2, 2))
  expect_equal(dim(mnix_res$staninput$transformed$Sigma_0_data), c(2, 2, 2))
})

context("check compatibility of new_ideal_adaptor_staninput with stanprograms")

test_that("new_ideal_adaptor_staninput accepts optional parameters via fix_parameters", {
  data <- make_minimal_staninput_data(cues = "cue1")

  res <- new_ideal_adaptor_staninput(
    exposure = data$exposure,
    test = data$test,
    cues = "cue1",
    category = "category",
    response = "response",
    group = "group",
    control = control_staninput(transform_type = "identity"),
    stanmodel = "NIX_ideal_adaptor",
    fix_parameters = list(
      mu_0 = list(0.1, 0.1),
      Sigma_0 = list(0.1, 0.1),
      lapse_rate = 0.05
    )
  )

  expect_equal(res$staninput$transformed$lapse_rate_data, 0.05)
  expect_equal(length(res$staninput$transformed$mu_0_data), 2)
  expect_equal(length(res$staninput$transformed$Sigma_0_data), 2)
})

test_that("NIW constructor supplies uncentered ss-style exposure inputs", {
  exposure <- data.frame(
    category = factor(c("A", "A"), levels = "A"),
    group = factor(c("g1", "g1"), levels = "g1"),
    response = factor(c("A", "A"), levels = "A"),
    cue1 = c(1, 2),
    cue2 = c(3, 4)
  )
  test <- data.frame(
    group = factor("g1", levels = "g1"),
    response = factor("A", levels = "A"),
    cue1 = 0.5,
    cue2 = 1.5
  )

  res <- new_ideal_adaptor_staninput(
    exposure = exposure,
    test = test,
    cues = c("cue1", "cue2"),
    category = "category",
    response = "response",
    group = "group",
    control = control_staninput(transform_type = "identity"),
    stanmodel = "NIW_ideal_adaptor"
  )

  expected <- crossprod(as.matrix(exposure[, c("cue1", "cue2")]))
  expect_equal(unname(res$staninput$transformed$x_ss_exposure[1, 1, , ]), unname(expected))
})

test_that("MNIX constructor supplies ss-style exposure inputs instead of covariance matrices", {
  data <- make_minimal_staninput_data(cues = c("cue1", "cue2"))

  res <- new_ideal_adaptor_staninput(
    exposure = data$exposure,
    test = data$test,
    cues = c("cue1", "cue2"),
    category = "category",
    response = "response",
    group = "group",
    control = control_staninput(transform_type = "identity"),
    stanmodel = "MNIX_ideal_adaptor"
  )

  expect_true("x_ss_exposure" %in% names(res$staninput$transformed))
  expect_false("x_cov_exposure" %in% names(res$staninput$transformed))
  expect_equal(dim(res$staninput$transformed$x_ss_exposure), c(2, 2, 2))
})

test_that("new_ideal_adaptor_staninput errors for empty test data", {
  data <- make_minimal_staninput_data(
    cues = "cue1",
    n_obs_exposure = 2L,
    n_obs_test = 0L,
    n_group = 1L
  )

  expect_error(
    new_ideal_adaptor_staninput(
      exposure = data$exposure,
      test = data$test,
      cues = "cue1",
      category = "category",
      response = "response",
      group = "group",
      control = control_staninput(transform_type = "identity"),
      stanmodel = "NIW_ideal_adaptor"
    ),
    "requires non-empty test data"
  )
})

test_that("new_ideal_adaptor_staninput accepts one or more test observations for a group", {
  for (n_obs_test in c(1L, 2L)) {
    data <- make_minimal_staninput_data(
      cues = "cue1",
      n_obs_exposure = 2L,
      n_obs_test = n_obs_test,
      n_group = 1L
    )

    res <- expect_no_error(
      new_ideal_adaptor_staninput(
        exposure = data$exposure,
        test = data$test,
        cues = "cue1",
        category = "category",
        response = "response",
        group = "group",
        control = control_staninput(transform_type = "identity"),
        stanmodel = "NIW_ideal_adaptor"
      )
    )

    expect_equal(res$staninput$transformed$N_test, n_obs_test)
  }
})

test_that("new_ideal_adaptor_staninput handles 0, 1, and multiple exposure observations for NIX", {
  skip_if_not_installed("rstan")

  for (n_obs in c(0L, 1L, 2L)) {
    data <- make_minimal_staninput_data(
      cues = "cue1",
      n_obs_exposure = n_obs,
      n_obs_test = 1L,
      n_group = 1L,
      n_category = 1L
    )
    res <- new_ideal_adaptor_staninput(
      exposure = data$exposure,
      test = data$test,
      cues = "cue1",
      category = "category",
      response = "response",
      group = "group",
      control = control_staninput(transform_type = "identity"),
      stanmodel = "NIX_ideal_adaptor"
    )

    expect_equal(res$staninput$transformed$N_exposure[1, 1], n_obs)
    expect_equal(dim(res$staninput$transformed$N_exposure), c(1, 1))

    fit <- rstan::sampling(
      object = MVBeliefUpdatr:::stanmodels[["NIX_ideal_adaptor"]],
      data = res$staninput$transformed,
      chains = 1,
      iter = 1,
      warmup = 0,
      refresh = 0,
      init = 0,
      algorithm = "Fixed_param"
    )

    expect_true(inherits(fit, "stanfit"))
  }
})

test_that("new_ideal_adaptor_staninput handles 0, 1, and multiple exposure observations for NIW", {
  skip_if_not_installed("rstan")

  for (n_obs in c(0L, 1L, 2L)) {
    data <- make_minimal_staninput_data(
      cues = c("cue1"),
      n_obs_exposure = n_obs,
      n_obs_test = 1L,
      n_group = 1L,
      n_category = 1L
    )
    res <- new_ideal_adaptor_staninput(
      exposure = data$exposure,
      test = data$test,
      cues = "cue1",
      category = "category",
      response = "response",
      group = "group",
      control = control_staninput(transform_type = "identity"),
      stanmodel = "NIW_ideal_adaptor"
    )

    expect_equal(res$staninput$transformed$N_exposure[1, 1], n_obs)
    expect_equal(dim(res$staninput$transformed$N_exposure), c(1, 1))

    fit <- rstan::sampling(
      object = MVBeliefUpdatr:::stanmodels[["NIW_ideal_adaptor"]],
      data = res$staninput$transformed,
      chains = 1,
      iter = 1,
      warmup = 0,
      refresh = 0,
      init = 0,
      algorithm = "Fixed_param"
    )

    expect_true(inherits(fit, "stanfit"))
  }
})

test_that("new_ideal_adaptor_staninput handles 0, 1, and multiple exposure observations for MNIX", {
  skip_if_not_installed("rstan")

  for (n_obs in c(0L, 1L, 2L)) {
    data <- make_minimal_staninput_data(
      cues = c("cue1", "cue2"),
      n_obs_exposure = n_obs,
      n_obs_test = 1L,
      n_group = 1L,
      n_category = 1L
    )
    res <- new_ideal_adaptor_staninput(
      exposure = data$exposure,
      test = data$test,
      cues = c("cue1", "cue2"),
      category = "category",
      response = "response",
      group = "group",
      control = control_staninput(transform_type = "identity"),
      stanmodel = "MNIX_ideal_adaptor"
    )

    expect_equal(res$staninput$transformed$N_exposure[1, 1], n_obs)
    expect_equal(dim(res$staninput$transformed$N_exposure), c(1, 1))

    fit <- rstan::sampling(
      object = MVBeliefUpdatr:::stanmodels[["MNIX_ideal_adaptor"]],
      data = res$staninput$transformed,
      chains = 1,
      iter = 1,
      warmup = 0,
      refresh = 0,
      init = 0,
      algorithm = "Fixed_param"
    )

    expect_true(inherits(fit, "stanfit"))
  }
})

test_that("new_ideal_adaptor_staninput produces inputs compatible with the NIW stan program", {
  skip_if_not_installed("rstan")

  data <- make_minimal_staninput_data(cues = "cue1")
  res <- new_ideal_adaptor_staninput(
    exposure = data$exposure,
    test = data$test,
    cues = "cue1",
    category = "category",
    response = "response",
    group = "group",
    control = control_staninput(transform_type = "identity"),
    stanmodel = "NIW_ideal_adaptor"
  )

  fit <- rstan::sampling(
    object = MVBeliefUpdatr:::stanmodels[["NIW_ideal_adaptor"]],
    data = res$staninput$transformed,
    chains = 1,
    iter = 1,
    warmup = 0,
    refresh = 0,
    init = 0,
    algorithm = "Fixed_param"
  )

  expect_true(inherits(fit, "stanfit"))
})

test_that("new_ideal_adaptor_staninput produces inputs compatible with the MNIX stan program", {
  skip_if_not_installed("rstan")

  data <- make_minimal_staninput_data(cues = c("cue1", "cue2"))
  res <- new_ideal_adaptor_staninput(
    exposure = data$exposure,
    test = data$test,
    cues = c("cue1", "cue2"),
    category = "category",
    response = "response",
    group = "group",
    control = control_staninput(transform_type = "identity"),
    stanmodel = "MNIX_ideal_adaptor"
  )

  fit <- rstan::sampling(
    object = MVBeliefUpdatr:::stanmodels[["MNIX_ideal_adaptor"]],
    data = res$staninput$transformed,
    chains = 1,
    iter = 1,
    warmup = 0,
    refresh = 0,
    init = 0,
    algorithm = "Fixed_param"
  )

  expect_true(inherits(fit, "stanfit"))
})

test_that("legacy constructor wrappers remain compatible with the new path", {
  data <- make_minimal_staninput_data(cues = "cue1")

  nix_input <- make_ideal_adaptor_stanfit_input(
    exposure = data$exposure,
    test = data$test,
    cues = "cue1",
    category = "category",
    response = "response",
    group = "group",
    control = control_staninput(transform_type = "identity"),
    stanmodel = "NIX_ideal_adaptor"
  )
  expect_true(is.ideal_adaptor_stanfit_input(nix_input))
  expect_true(is.data.frame(nix_input$data))
  expect_true(is.list(nix_input$staninput$transformed))
  expect_true(is.list(nix_input$staninput$untransformed))

  niw_input <- make_ideal_adaptor_stanfit_input(
    exposure = data$exposure,
    test = data$test,
    cues = "cue1",
    category = "category",
    response = "response",
    group = "group",
    control = control_staninput(transform_type = "standardize"),
    stanmodel = "NIW_ideal_adaptor"
  )
  expect_true(is.ideal_adaptor_stanfit_input(niw_input))
  expect_true("x_mean_exposure" %in% names(niw_input$staninput$untransformed))

  multi_cue_data <- make_minimal_staninput_data(cues = c("cue1", "cue2"))
  mnix_input <- make_ideal_adaptor_stanfit_input(
    exposure = multi_cue_data$exposure,
    test = multi_cue_data$test,
    cues = c("cue1", "cue2"),
    category = "category",
    response = "response",
    group = "group",
    control = control_staninput(transform_type = "identity"),
    stanmodel = "MNIX_ideal_adaptor"
  )
  expect_true(is.ideal_adaptor_stanfit_input(mnix_input))
})

test_that("invalid transform types fail early in legacy wrappers", {
  data <- make_minimal_staninput_data(cues = "cue1")

  expect_error(
    make_ideal_adaptor_stanfit_input(
      exposure = data$exposure,
      test = data$test,
      cues = "cue1",
      category = "category",
      response = "response",
      group = "group",
      control = control_staninput(transform_type = "other"),
      stanmodel = "NIW_ideal_adaptor"
    ),
    "transform_type"
  )
})

test_that("model-specific cue requirements are enforced by legacy wrappers", {
  data <- make_minimal_staninput_data(cues = c("cue1", "cue2"))

  expect_error(
    make_ideal_adaptor_stanfit_input(
      exposure = data$exposure,
      test = data$test,
      cues = c("cue1", "cue2"),
      category = "category",
      response = "response",
      group = "group",
      control = control_staninput(transform_type = "identity"),
      stanmodel = "NIX_ideal_adaptor"
    ),
    "requires exactly one cue"
  )

  expect_error(
    make_ideal_adaptor_stanfit_input(
      exposure = data$exposure,
      test = data$test,
      cues = "cue1",
      category = "category",
      response = "response",
      group = "group",
      control = control_staninput(transform_type = "identity"),
      stanmodel = "MNIX_ideal_adaptor"
    ),
    "requires at least two cues"
  )
})

test_that("invalid lapse rates are rejected by legacy wrappers", {
  data <- make_minimal_staninput_data(cues = "cue1")

  expect_error(
    make_ideal_adaptor_stanfit_input(
      exposure = data$exposure,
      test = data$test,
      cues = "cue1",
      category = "category",
      response = "response",
      group = "group",
      lapse_rate = 1.5,
      control = control_staninput(transform_type = "identity"),
      stanmodel = "NIW_ideal_adaptor"
    ),
    "between 0 and 1"
  )
})
