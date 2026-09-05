test_that("new_ideal_adaptor_stanfit_input rejects missing or empty cues", {
  data <- build_shifted_prior_fit_example(
    cues = c("VOT", "f0_semitones"),
    n_exposure_per_condition = 20L,
    n_test_per_condition = 10L,
    seed = 321L
  )

  expect_error(
    do.call(
      new_ideal_adaptor_stanfit_input,
      list(
        exposure = data$exposure,
        test = data$test,
        category = "category",
        response = "response",
        group = "group",
        group.unique = "Condition",
        control = control_staninput(transform_type = "identity"),
        stanmodel = "MNIX_ideal_adaptor"
      )
    ),
    regexp = "Expected x to be a non-NA character",
    fixed = TRUE
  )

  expect_error(
    new_ideal_adaptor_stanfit_input(
      exposure = data$exposure,
      test = data$test,
      cues = character(0),
      category = "category",
      response = "response",
      group = "group",
      group.unique = "Condition",
      control = control_staninput(transform_type = "identity"),
      stanmodel = "MNIX_ideal_adaptor"
    ),
    regexp = "Expected x to be a non-NA character",
    fixed = TRUE
  )
})

test_that("new_ideal_adaptor_stanfit_input returns a structured fit-input object", {
  exposure <- data.frame(
    category = factor(c("A", "B", "A", "B")),
    group = factor(c("g1", "g1", "g2", "g2")),
    cue1 = c(0.1, 0.4, 0.2, 0.5),
    cue2 = c(0.2, 0.3, 0.1, 0.6)
  )
  test <- data.frame(
    response = factor(c("A", "B", "A", "B")),
    group = factor(c("g1", "g2", "g1", "g2")),
    cue1 = c(0.15, 0.45, 0.25, 0.55),
    cue2 = c(0.25, 0.35, 0.15, 0.65)
  )

  input <- new_ideal_adaptor_stanfit_input(
    exposure = exposure,
    test = test,
    cues = c("cue1", "cue2"),
    category = "category",
    response = "response",
    group = "group",
    stanmodel = "NIW_ideal_adaptor"
  )

  expect_true(S7::S7_inherits(input, IdealAdaptorStanfitInput))
  expect_s3_class(input@data, "data.frame")
  expect_true(S7::S7_inherits(input@staninput, IdealAdaptorStaninput))
  expect_true(S7::S7_inherits(input@transform_information, MVBU_TransformInformation))
})

test_that("new_ideal_adaptor_stanfit_input preserves grouped exposure data", {
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

  res <- new_ideal_adaptor_stanfit_input(
    exposure = exposure,
    test = test,
    cues = "cue1",
    category = "category",
    response = "response",
    group = "group",
    group.unique = "Condition",
    control = control_staninput(transform_type = "identity"),
    stanmodel = "NIW_ideal_adaptor"
  )

  expect_true("Condition" %in% names(res@data))
  expect_true(all(res@data$Condition[res@data$Phase == "exposure"] == "baseline"))
})

test_that("new_ideal_adaptor_stanfit_input handles more than two categories in the input", {
  exposure <- data.frame(
    category = factor(c("A", "B", "C", "A", "B", "C"), levels = c("A", "B", "C")),
    group = factor(c("g1", "g1", "g1", "g2", "g2", "g2"), levels = c("g1", "g2")),
    cue1 = c(0.1, 0.4, 0.7, 0.2, 0.5, 0.8),
    cue2 = c(0.2, 0.3, 0.6, 0.1, 0.4, 0.9)
  )
  test <- data.frame(
    response = factor(c("A", "B", "C"), levels = c("A", "B", "C")),
    group = factor(c("g1", "g2", "g1"), levels = c("g1", "g2")),
    cue1 = c(0.15, 0.55, 0.75),
    cue2 = c(0.25, 0.45, 0.65)
  )

  res <- new_ideal_adaptor_stanfit_input(
    exposure = exposure,
    test = test,
    cues = c("cue1", "cue2"),
    category = "category",
    response = "response",
    group = "group",
    control = control_staninput(transform_type = "identity"),
    stanmodel = "NIW_ideal_adaptor"
  )

  expect_true(S7::S7_inherits(res, IdealAdaptorStanfitInput))
  expect_true(S7::S7_inherits(res@staninput, IdealAdaptorStaninput))
  expect_equal(as.matrix(res@staninput@values$N_exposure), matrix(1L, nrow = 3, ncol = 2))
  expect_equal(res@staninput@values$N_test, 3L)
})

test_that("new_ideal_adaptor_stanfit_input rejects empty test data", {
  data <- make_minimal_staninput_data(
    cues = "cue1",
    n_obs_exposure = 2L,
    n_obs_test = 0L,
    n_group = 1L
  )

  expect_error(
    new_ideal_adaptor_stanfit_input(
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

test_that("fixed_parameters are validated against the category labels", {
  exposure <- data.frame(
    category = factor(c("A", "B")),
    group = factor(c("g1", "g1")),
    cue1 = c(0, 1),
    cue2 = c(0, 1)
  )
  test <- data.frame(
    response = factor(c("A", "B")),
    group = factor(c("g1", "g1")),
    cue1 = c(0.1, 0.9),
    cue2 = c(0.1, 0.9)
  )

  expect_error(
    new_ideal_adaptor_stanfit_input(
      exposure = exposure,
      test = test,
      cues = c("cue1", "cue2"),
      category = "category",
      response = "response",
      group = "group",
      fixed_parameters = list(
        mu_0 = list(B = c(0, 0)),
        Sigma_0 = list(B = diag(2))
      ),
      control = control_staninput(transform_type = "identity"),
      stanmodel = "NIW_ideal_adaptor"
    ),
    "named list with names matching the categories"
  )
})

test_that("fixed_parameters support single, pairwise, and full combinations across Stan families", {
  one_cue_data <- make_minimal_staninput_data(cues = "cue1")
  two_cue_data <- make_minimal_staninput_data(cues = c("cue1", "cue2"))

  nix_fixed_sets <- list(
    list(lapse_rate = 0.1),
    list(mu_0 = list("A1" = 0.1, "A2" = 0.2)),
    list(Sigma_0 = list("A1" = matrix(0.1, nrow = 1, ncol = 1), "A2" = matrix(0.2, nrow = 1, ncol = 1))),
    list(lapse_rate = 0.1, mu_0 = list("A1" = 0.1, "A2" = 0.2)),
    list(lapse_rate = 0.1, Sigma_0 = list("A1" = matrix(0.1, nrow = 1, ncol = 1), "A2" = matrix(0.2, nrow = 1, ncol = 1))),
    list(mu_0 = list("A1" = 0.1, "A2" = 0.2), Sigma_0 = list("A1" = matrix(0.1, nrow = 1, ncol = 1), "A2" = matrix(0.2, nrow = 1, ncol = 1))),
    list(lapse_rate = 0.1, mu_0 = list("A1" = 0.1, "A2" = 0.2), Sigma_0 = list("A1" = matrix(0.1, nrow = 1, ncol = 1), "A2" = matrix(0.2, nrow = 1, ncol = 1)))
  )

  mnix_niw_fixed_sets <- list(
    list(lapse_rate = 0.1),
    list(mu_0 = list("A1" = c(0.1, 0.2), "A2" = c(0.2, 0.3))),
    list(Sigma_0 = list("A1" = diag(2), "A2" = diag(2))),
    list(lapse_rate = 0.1, mu_0 = list("A1" = c(0.1, 0.2), "A2" = c(0.2, 0.3))),
    list(lapse_rate = 0.1, Sigma_0 = list("A1" = diag(2), "A2" = diag(2))),
    list(mu_0 = list("A1" = c(0.1, 0.2), "A2" = c(0.2, 0.3)), Sigma_0 = list("A1" = diag(2), "A2" = diag(2))),
    list(lapse_rate = 0.1, mu_0 = list("A1" = c(0.1, 0.2), "A2" = c(0.2, 0.3)), Sigma_0 = list("A1" = diag(2), "A2" = diag(2)))
  )

  for (fixed_parameters in nix_fixed_sets) {
    nix_res <- new_ideal_adaptor_stanfit_input(
      exposure = one_cue_data$exposure,
      test = one_cue_data$test,
      cues = "cue1",
      category = "category",
      response = "response",
      group = "group",
      control = control_staninput(transform_type = "identity"),
      stanmodel = "NIX_ideal_adaptor",
      fixed_parameters = fixed_parameters
    )

    expect_true(S7::S7_inherits(nix_res, IdealAdaptorStanfitInput))
    expect_equal(nix_res@staninput@values$lapse_rate_known, if (is.null(fixed_parameters$lapse_rate)) 0 else 1)
    expect_equal(nix_res@staninput@values$mu_0_known, if (is.null(fixed_parameters$mu_0)) 0 else 1)
    expect_equal(nix_res@staninput@values$Sigma_0_known, if (is.null(fixed_parameters$Sigma_0)) 0 else 1)
  }

  for (fixed_parameters in mnix_niw_fixed_sets) {
    mnix_res <- new_ideal_adaptor_stanfit_input(
      exposure = two_cue_data$exposure,
      test = two_cue_data$test,
      cues = c("cue1", "cue2"),
      category = "category",
      response = "response",
      group = "group",
      control = control_staninput(transform_type = "identity"),
      stanmodel = "MNIX_ideal_adaptor",
      fixed_parameters = fixed_parameters
    )

    expect_true(S7::S7_inherits(mnix_res, IdealAdaptorStanfitInput))
    expect_equal(mnix_res@staninput@values$lapse_rate_known, if (is.null(fixed_parameters$lapse_rate)) 0 else 1)
    expect_equal(mnix_res@staninput@values$mu_0_known, if (is.null(fixed_parameters$mu_0)) 0 else 1)
    expect_equal(mnix_res@staninput@values$Sigma_0_known, if (is.null(fixed_parameters$Sigma_0)) 0 else 1)

    niw_res <- new_ideal_adaptor_stanfit_input(
      exposure = two_cue_data$exposure,
      test = two_cue_data$test,
      cues = c("cue1", "cue2"),
      category = "category",
      response = "response",
      group = "group",
      control = control_staninput(transform_type = "identity"),
      stanmodel = "NIW_ideal_adaptor",
      fixed_parameters = fixed_parameters
    )

    expect_true(S7::S7_inherits(niw_res, IdealAdaptorStanfitInput))
    expect_equal(niw_res@staninput@values$lapse_rate_known, if (is.null(fixed_parameters$lapse_rate)) 0 else 1)
    expect_equal(niw_res@staninput@values$mu_0_known, if (is.null(fixed_parameters$mu_0)) 0 else 1)
    expect_equal(niw_res@staninput@values$Sigma_0_known, if (is.null(fixed_parameters$Sigma_0)) 0 else 1)
  }
})

test_that("fixed_parameters reject invalid dimensionality or incomplete information", {
  two_cue_data <- make_minimal_staninput_data(cues = c("cue1", "cue2"))

  expect_error(
    new_ideal_adaptor_stanfit_input(
      exposure = two_cue_data$exposure,
      test = two_cue_data$test,
      cues = c("cue1", "cue2"),
      category = "category",
      response = "response",
      group = "group",
      control = control_staninput(transform_type = "identity"),
      stanmodel = "NIW_ideal_adaptor",
      fixed_parameters = list(mu_0 = list("A1" = 0.1, "A2" = c(0.1, 0.2, 0.3)))
    ),
    "mu_0 must be a named list with each entry a numeric vector of length 2"
  )

  expect_error(
    new_ideal_adaptor_stanfit_input(
      exposure = two_cue_data$exposure,
      test = two_cue_data$test,
      cues = c("cue1", "cue2"),
      category = "category",
      response = "response",
      group = "group",
      control = control_staninput(transform_type = "identity"),
      stanmodel = "NIW_ideal_adaptor",
      fixed_parameters = list(Sigma_0 = list("A1" = diag(1), "A2" = diag(2)))
    ),
    "Sigma_0 must be a named list with each entry a square matrix of dimension 2 x 2"
  )

  expect_error(
    new_ideal_adaptor_stanfit_input(
      exposure = two_cue_data$exposure,
      test = two_cue_data$test,
      cues = c("cue1", "cue2"),
      category = "category",
      response = "response",
      group = "group",
      control = control_staninput(transform_type = "identity"),
      stanmodel = "MNIX_ideal_adaptor",
      fixed_parameters = list(mu_0 = list("A1" = c(0.1, 0.2), "A2" = c(0.1, 0.2), "A3" = c(0.1, 0.2)))
    ),
    "mu_0 must be a named list with names matching the categories present in the exposure and test data"
  )
})

test_that("stanmodel compatibility checks enforce the matching Stan program", {
  for (n_obs_exposure in c(0L, 1L, 3L)) {
    data <- make_minimal_staninput_data(
      cues = "cue1",
      n_obs_exposure = n_obs_exposure,
      n_obs_test = 2L,
      n_group = 1L,
      n_category = 2L
    )

    nix_input <- new_ideal_adaptor_stanfit_input(
      exposure = data$exposure,
      test = data$test,
      cues = "cue1",
      category = "category",
      response = "response",
      group = "group",
      control = control_staninput(transform_type = "identity"),
      stanmodel = "NIX_ideal_adaptor"
    )
    expect_staninput_structure(
      nix_input,
      expected_class = NIX_IdealAdaptorStaninput,
      required_names = c("N_exposure", "N_test", "x_mean_exposure", "x_sd_exposure"),
      forbidden_names = c("x_ss_exposure")
    )

    data_two_cues <- make_minimal_staninput_data(
      cues = c("cue1", "cue2"),
      n_obs_exposure = n_obs_exposure,
      n_obs_test = 2L,
      n_group = 1L,
      n_category = 2L
    )

    expect_error(
      new_ideal_adaptor_stanfit_input(
        exposure = data_two_cues$exposure,
        test = data_two_cues$test,
        cues = c("cue1", "cue2"),
        category = "category",
        response = "response",
        group = "group",
        control = control_staninput(transform_type = "identity"),
        stanmodel = "NIX_ideal_adaptor"
      ),
      "requires exactly one cue"
    )
  }

  for (n_obs_exposure in c(0L, 1L, 3L)) {
    data <- make_minimal_staninput_data(
      cues = c("cue1", "cue2"),
      n_obs_exposure = n_obs_exposure,
      n_obs_test = 2L,
      n_group = 1L,
      n_category = 2L
    )

    mnix_input <- new_ideal_adaptor_stanfit_input(
      exposure = data$exposure,
      test = data$test,
      cues = c("cue1", "cue2"),
      category = "category",
      response = "response",
      group = "group",
      control = control_staninput(transform_type = "identity"),
      stanmodel = "MNIX_ideal_adaptor"
    )
    expect_staninput_structure(
      mnix_input,
      expected_class = MNIX_IdealAdaptorStaninput,
      required_names = c("N_exposure", "N_test", "x_mean_exposure", "x_ss_exposure", "p_cat"),
      forbidden_names = c("x_sd_exposure")
    )

    data_one_cue <- make_minimal_staninput_data(
      cues = "cue1",
      n_obs_exposure = n_obs_exposure,
      n_obs_test = 2L,
      n_group = 1L,
      n_category = 2L
    )

    expect_error(
      new_ideal_adaptor_stanfit_input(
        exposure = data_one_cue$exposure,
        test = data_one_cue$test,
        cues = "cue1",
        category = "category",
        response = "response",
        group = "group",
        control = control_staninput(transform_type = "identity"),
        stanmodel = "MNIX_ideal_adaptor"
      ),
      "requires at least two cues"
    )
  }

  for (n_obs_exposure in c(0L, 1L, 3L)) {
    for (cues in list("cue1", c("cue1", "cue2"), c("cue1", "cue2", "cue3"))) {
      data <- make_minimal_staninput_data(
        cues = cues,
        n_obs_exposure = n_obs_exposure,
        n_obs_test = 2L,
        n_group = 1L,
        n_category = 2L
      )

      niw_input <- new_ideal_adaptor_stanfit_input(
        exposure = data$exposure,
        test = data$test,
        cues = cues,
        category = "category",
        response = "response",
        group = "group",
        control = control_staninput(transform_type = "identity"),
        stanmodel = "NIW_ideal_adaptor"
      )
      expect_staninput_structure(
        niw_input,
        expected_class = NIW_IdealAdaptorStaninput,
        required_names = c("N_exposure", "N_test", "x_mean_exposure", "x_ss_exposure"),
        forbidden_names = c("x_sd_exposure")
      )
    }
  }
})

test_that("NIX supports one cue, rejects 2 or more cues, and matches the NIX Stan input structure", {
  for (n_obs_exposure in c(0L, 1L, 3L)) {
    data <- make_minimal_staninput_data(
      cues = "cue1",
      n_obs_exposure = n_obs_exposure,
      n_obs_test = 2L,
      n_group = 1L,
      n_category = 2L
    )

    input <- new_ideal_adaptor_stanfit_input(
      exposure = data$exposure,
      test = data$test,
      cues = "cue1",
      category = "category",
      response = "response",
      group = "group",
      control = control_staninput(transform_type = "identity"),
      stanmodel = "NIX_ideal_adaptor"
    )

    expected_counts <- tabulate(rep(seq_len(2), length.out = n_obs_exposure), nbins = 2)

    expect_staninput_structure(
      input,
      expected_class = NIX_IdealAdaptorStaninput,
      required_names = c("N_exposure", "N_test", "x_mean_exposure", "x_sd_exposure"),
      forbidden_names = c("x_ss_exposure")
    )
    expect_equal(as.vector(input@staninput@values$N_exposure), expected_counts)
    expect_equal(input@staninput@values$N_test, 2L)
  }

  data <- make_minimal_staninput_data(cues = c("cue1", "cue2"), n_obs_exposure = 1L, n_obs_test = 2L)
  expect_error(
    new_ideal_adaptor_stanfit_input(
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
})

test_that("MNIX supports 2 or more cues, rejects 1 cue, and matches the MNIX Stan input structure", {
  for (n_obs_exposure in c(0L, 1L, 3L)) {
    data <- make_minimal_staninput_data(
      cues = c("cue1", "cue2"),
      n_obs_exposure = n_obs_exposure,
      n_obs_test = 2L,
      n_group = 1L,
      n_category = 2L
    )

    input <- new_ideal_adaptor_stanfit_input(
      exposure = data$exposure,
      test = data$test,
      cues = c("cue1", "cue2"),
      category = "category",
      response = "response",
      group = "group",
      control = control_staninput(transform_type = "identity"),
      stanmodel = "MNIX_ideal_adaptor"
    )

    expected_counts <- tabulate(rep(seq_len(2), length.out = n_obs_exposure), nbins = 2)

    expect_staninput_structure(
      input,
      expected_class = MNIX_IdealAdaptorStaninput,
      required_names = c("N_exposure", "N_test", "x_mean_exposure", "x_ss_exposure", "p_cat"),
      forbidden_names = c("x_sd_exposure")
    )
    expect_equal(as.vector(input@staninput@values$N_exposure), expected_counts)
    expect_equal(input@staninput@values$N_test, 2L)
  }

  data <- make_minimal_staninput_data(cues = "cue1", n_obs_exposure = 1L, n_obs_test = 2L)
  expect_error(
    new_ideal_adaptor_stanfit_input(
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

test_that("NIW supports 1, 2, and 2+ cues and matches the NIW Stan input structure", {
  for (cues in list("cue1", c("cue1", "cue2"), c("cue1", "cue2", "cue3"))) {
    for (n_obs_exposure in c(0L, 1L, 3L)) {
      data <- make_minimal_staninput_data(
        cues = cues,
        n_obs_exposure = n_obs_exposure,
        n_obs_test = 2L,
        n_group = 1L,
        n_category = 2L
      )

      input <- new_ideal_adaptor_stanfit_input(
        exposure = data$exposure,
        test = data$test,
        cues = cues,
        category = "category",
        response = "response",
        group = "group",
        control = control_staninput(transform_type = "identity"),
        stanmodel = "NIW_ideal_adaptor"
      )

      expected_counts <- tabulate(rep(seq_len(2), length.out = n_obs_exposure), nbins = 2)

      expect_staninput_structure(
        input,
        expected_class = NIW_IdealAdaptorStaninput,
        required_names = c("N_exposure", "N_test", "x_mean_exposure", "x_ss_exposure"),
        forbidden_names = c("x_sd_exposure")
      )
      expect_equal(as.vector(input@staninput@values$N_exposure), expected_counts)
      expect_equal(input@staninput@values$N_test, 2L)
    }
  }
})

test_that("group.unique correctly checks identity and simplifies staninput", {
  # Case 1: Identical exposure stats across groups within unique groups
  exp1 <- data.frame(
    Condition = factor(c("c1", "c1", "c1", "c1", "c2", "c2")),
    group = factor(c("g1", "g1", "g2", "g2", "g3", "g3")),
    category = factor(c("A", "B", "A", "B", "A", "B")),
    cue1 = c(1, 2, 1, 2, 5, 6),
    cue2 = c(3, 4, 3, 4, 7, 8)
  )
  test1 <- data.frame(
    Condition = factor(c("c1", "c1", "c2")),
    group = factor(c("g1", "g2", "g3")),
    response = factor(c("A", "B", "A")),
    cue1 = c(1.1, 2.1, 5.1),
    cue2 = c(3.1, 4.1, 7.1)
  )

  input1 <- new_ideal_adaptor_stanfit_input(
    exposure = exp1,
    test = test1,
    cues = c("cue1", "cue2"),
    category = "category",
    response = "response",
    group = "group",
    group.unique = "Condition",
    check_unique_group_identity = TRUE,
    control = control_staninput(transform_type = "identity"),
    stanmodel = "NIW_ideal_adaptor"
  )

  expect_equal(input1@staninput@values$L, 2L)
  expect_equal(dim(input1@staninput@values$N_exposure), c(2L, 2L))
  expect_equal(dim(input1@staninput@values$x_mean_exposure), c(2L, 2L, 2L))
  expect_equal(dim(input1@staninput@values$x_ss_exposure), c(2L, 2L, 2L, 2L))
  expect_equal(
    get_labels(input1)$group,
    c("c1", "c2")
  )
  # Test rows map to c1 (index 1) and c2 (index 2)
  expect_equal(input1@staninput@values$y_test, c(1L, 1L, 2L))

  # Case 2: Mismatching exposure stats with check_unique_group_identity = TRUE
  exp_mismatch <- exp1
  exp_mismatch$cue1[3] <- 99 # g2 in c1 now has different mean

  expect_error(
    new_ideal_adaptor_stanfit_input(
      exposure = exp_mismatch,
      test = test1,
      cues = c("cue1", "cue2"),
      category = "category",
      response = "response",
      group = "group",
      group.unique = "Condition",
      check_unique_group_identity = TRUE,
      control = control_staninput(transform_type = "identity"),
      stanmodel = "NIW_ideal_adaptor"
    ),
    regexp = "Non-identical exposure sufficient statistics found within unique group(s): c1",
    fixed = TRUE
  )

  # Case 3: Mismatching exposure stats with check_unique_group_identity = FALSE
  input_no_check <- new_ideal_adaptor_stanfit_input(
    exposure = exp_mismatch,
    test = test1,
    cues = c("cue1", "cue2"),
    category = "category",
    response = "response",
    group = "group",
    group.unique = "Condition",
    check_unique_group_identity = FALSE,
    control = control_staninput(transform_type = "identity"),
    stanmodel = "NIW_ideal_adaptor"
  )
  expect_equal(input_no_check@staninput@values$L, 2L)
  expect_equal(get_labels(input_no_check)$group, c("c1", "c2"))
})


