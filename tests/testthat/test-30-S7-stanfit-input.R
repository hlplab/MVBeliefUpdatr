test_that("new_ideal_adaptor_stanfit_input rejects missing or empty cues", {
  data <- example_exposure_test_data(
    model_family = "MNIX",
    n_cues = 2L,
    seed = 321L
  )

  expect_error(
    do.call(
      new_ideal_adaptor_stanfit_input,
      list(
        exposure = data$exposure,
        test = data$test,
        category = "category",
        response = "response_category",
        group = "group",
        group_unique = "Condition",
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
      response = "response_category",
      group = "group",
      group_unique = "Condition",
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
    group_unique = "Condition",
    control = control_staninput(transform_type = "identity"),
    stanmodel = "NIW_ideal_adaptor"
  )

  expect_true("group_unique" %in% names(res@data))
  expect_true(all(res@data$group_unique[res@data$Phase == "exposure"] == "baseline"))
  orig_data <- get_data(res, original_names = TRUE)
  expect_true("Condition" %in% names(orig_data))
  expect_true(all(orig_data$Condition[orig_data$Phase == "exposure"] == "baseline"))
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
  data <- example_exposure_test_data(model_family = "NIW", n_cues = 1L)
  data$test <- data$test[0, , drop = FALSE]

  expect_error(
    new_ideal_adaptor_stanfit_input(
      exposure = data$exposure,
      test = data$test,
      cues = "VOT",
      category = "category",
      response = "response_category",
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
    response_category = factor(c("A", "B")),
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
      response = "response_category",
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

test_that(
  "fixed_parameters support single, pairwise, and full combinations",
  {
  one_cue_data <- example_exposure_test_data(
    model_family = "NIX",
    n_cues = 1L
  )
  two_cue_data <- example_exposure_test_data(
    model_family = "MNIX",
    n_cues = 2L
  )

  nix_fixed_sets <- list(
    list(lapse_rate = 0.1),
    list(mu_0 = list("/d/" = 0.1, "/t/" = 0.2)),
    list(
      Sigma_0 = list(
        "/d/" = matrix(0.1, nrow = 1, ncol = 1),
        "/t/" = matrix(0.2, nrow = 1, ncol = 1)
      )
    ),
    list(lapse_rate = 0.1, mu_0 = list("/d/" = 0.1, "/t/" = 0.2)),
    list(
      lapse_rate = 0.1,
      Sigma_0 = list(
        "/d/" = matrix(0.1, nrow = 1, ncol = 1),
        "/t/" = matrix(0.2, nrow = 1, ncol = 1)
      )
    ),
    list(
      mu_0 = list("/d/" = 0.1, "/t/" = 0.2),
      Sigma_0 = list(
        "/d/" = matrix(0.1, nrow = 1, ncol = 1),
        "/t/" = matrix(0.2, nrow = 1, ncol = 1)
      )
    ),
    list(
      lapse_rate = 0.1,
      mu_0 = list("/d/" = 0.1, "/t/" = 0.2),
      Sigma_0 = list(
        "/d/" = matrix(0.1, nrow = 1, ncol = 1),
        "/t/" = matrix(0.2, nrow = 1, ncol = 1)
      )
    )
  )

  mnix_niw_fixed_sets <- list(
    list(lapse_rate = 0.1),
    list(mu_0 = list("/d/" = c(0.1, 0.2), "/t/" = c(0.2, 0.3))),
    list(Sigma_0 = list("/d/" = diag(2), "/t/" = diag(2))),
    list(
      lapse_rate = 0.1,
      mu_0 = list("/d/" = c(0.1, 0.2), "/t/" = c(0.2, 0.3))
    ),
    list(lapse_rate = 0.1, Sigma_0 = list("/d/" = diag(2), "/t/" = diag(2))),
    list(
      mu_0 = list("/d/" = c(0.1, 0.2), "/t/" = c(0.2, 0.3)),
      Sigma_0 = list("/d/" = diag(2), "/t/" = diag(2))
    ),
    list(
      lapse_rate = 0.1,
      mu_0 = list("/d/" = c(0.1, 0.2), "/t/" = c(0.2, 0.3)),
      Sigma_0 = list("/d/" = diag(2), "/t/" = diag(2))
    )
  )

  for (fixed_parameters in nix_fixed_sets) {
    nix_res <- new_ideal_adaptor_stanfit_input(
      exposure = one_cue_data$exposure,
      test = one_cue_data$test,
      cues = "VOT",
      category = "category",
      response = "response_category",
      group = "group",
      control = control_staninput(transform_type = "identity"),
      stanmodel = "NIX_ideal_adaptor",
      fixed_parameters = fixed_parameters
    )

    expect_true(S7::S7_inherits(nix_res, IdealAdaptorStanfitInput))
    expect_equal(
      nix_res@staninput@values$lapse_rate_known,
      if (is.null(fixed_parameters$lapse_rate)) 0 else 1
    )
    expect_equal(
      nix_res@staninput@values$mu_0_known,
      if (is.null(fixed_parameters$mu_0)) 0 else 1
    )
    expect_equal(
      nix_res@staninput@values$Sigma_0_known,
      if (is.null(fixed_parameters$Sigma_0)) 0 else 1
    )
  }

  for (fixed_parameters in mnix_niw_fixed_sets) {
    mnix_res <- new_ideal_adaptor_stanfit_input(
      exposure = two_cue_data$exposure,
      test = two_cue_data$test,
      cues = c("VOT", "f0_semitones"),
      category = "category",
      response = "response_category",
      group = "group",
      control = control_staninput(transform_type = "identity"),
      stanmodel = "MNIX_ideal_adaptor",
      fixed_parameters = fixed_parameters
    )

    expect_true(S7::S7_inherits(mnix_res, IdealAdaptorStanfitInput))
    expect_equal(
      mnix_res@staninput@values$lapse_rate_known,
      if (is.null(fixed_parameters$lapse_rate)) 0 else 1
    )
    expect_equal(
      mnix_res@staninput@values$mu_0_known,
      if (is.null(fixed_parameters$mu_0)) 0 else 1
    )
    expect_equal(
      mnix_res@staninput@values$Sigma_0_known,
      if (is.null(fixed_parameters$Sigma_0)) 0 else 1
    )

    niw_res <- new_ideal_adaptor_stanfit_input(
      exposure = two_cue_data$exposure,
      test = two_cue_data$test,
      cues = c("VOT", "f0_semitones"),
      category = "category",
      response = "response_category",
      group = "group",
      control = control_staninput(transform_type = "identity"),
      stanmodel = "NIW_ideal_adaptor",
      fixed_parameters = fixed_parameters
    )

    expect_true(S7::S7_inherits(niw_res, IdealAdaptorStanfitInput))
    expect_equal(
      niw_res@staninput@values$lapse_rate_known,
      if (is.null(fixed_parameters$lapse_rate)) 0 else 1
    )
    expect_equal(
      niw_res@staninput@values$mu_0_known,
      if (is.null(fixed_parameters$mu_0)) 0 else 1
    )
    expect_equal(
      niw_res@staninput@values$Sigma_0_known,
      if (is.null(fixed_parameters$Sigma_0)) 0 else 1
    )
  }
})

test_that(
  "fixed_parameters reject invalid dimensionality or incomplete information",
  {
  two_cue_data <- example_exposure_test_data(
    model_family = "MNIX",
    n_cues = 2L
  )

  expect_error(
    new_ideal_adaptor_stanfit_input(
      exposure = two_cue_data$exposure,
      test = two_cue_data$test,
      cues = c("VOT", "f0_semitones"),
      category = "category",
      response = "response_category",
      group = "group",
      control = control_staninput(transform_type = "identity"),
      stanmodel = "NIW_ideal_adaptor",
      fixed_parameters = list(
        mu_0 = list("/d/" = 0.1, "/t/" = c(0.1, 0.2, 0.3))
      )
    ),
    "mu_0 must be a named list with each entry a numeric vector of length 2"
  )

  expect_error(
    new_ideal_adaptor_stanfit_input(
      exposure = two_cue_data$exposure,
      test = two_cue_data$test,
      cues = c("VOT", "f0_semitones"),
      category = "category",
      response = "response_category",
      group = "group",
      control = control_staninput(transform_type = "identity"),
      stanmodel = "NIW_ideal_adaptor",
      fixed_parameters = list(
        Sigma_0 = list("/d/" = diag(1), "/t/" = diag(2))
      )
    ),
    "Sigma_0 must be a named list with each entry a square matrix"
  )

  expect_error(
    new_ideal_adaptor_stanfit_input(
      exposure = two_cue_data$exposure,
      test = two_cue_data$test,
      cues = c("VOT", "f0_semitones"),
      category = "category",
      response = "response_category",
      group = "group",
      control = control_staninput(transform_type = "identity"),
      stanmodel = "MNIX_ideal_adaptor",
      fixed_parameters = list(
        mu_0 = list(
          "/d/" = c(0.1, 0.2),
          "/t/" = c(0.1, 0.2),
          "/b/" = c(0.1, 0.2)
        )
      )
    ),
    "mu_0 must be a named list with names matching the categories"
  )
})

test_that("stanmodel compatibility checks enforce the matching Stan program", {
  base_1cue <- example_exposure_test_data(model_family = "NIX", n_cues = 1L)
  base_2cue <- example_exposure_test_data(model_family = "MNIX", n_cues = 2L)
  base_3cue <- example_exposure_test_data(model_family = "NIW", n_cues = 3L)

  base_1cue$exposure$group <- "g1"
  base_1cue$test$group <- "g1"
  base_2cue$exposure$group <- "g1"
  base_2cue$test$group <- "g1"
  base_3cue$exposure$group <- "g1"
  base_3cue$test$group <- "g1"

  test_1cue <- base_1cue$test[c(1, 41), ]
  test_2cue <- base_2cue$test[c(1, 41), ]
  test_3cue <- base_3cue$test[c(1, 41), ]

  for (n_obs_exposure in c(0L, 1L, 3L)) {
    exp_1cue <- if (n_obs_exposure == 0L) {
      base_1cue$exposure[0L, ]
    } else {
      base_1cue$exposure[seq_len(n_obs_exposure), ]
    }
    exp_2cue <- if (n_obs_exposure == 0L) {
      base_2cue$exposure[0L, ]
    } else {
      base_2cue$exposure[seq_len(n_obs_exposure), ]
    }

    nix_input <- new_ideal_adaptor_stanfit_input(
      exposure = exp_1cue,
      test = test_1cue,
      cues = "VOT",
      category = "category",
      response = "response_category",
      group = "group",
      control = control_staninput(transform_type = "identity"),
      stanmodel = "NIX_ideal_adaptor"
    )
    expect_staninput_structure(
      nix_input,
      expected_class = NIX_IdealAdaptorStaninput,
      required_names = c(
        "N_exposure", "N_test", "x_mean_exposure", "x_ss_exposure"
      ),
      forbidden_names = c("x_sd_exposure")
    )

    expect_error(
      new_ideal_adaptor_stanfit_input(
        exposure = exp_2cue,
        test = test_2cue,
        cues = c("VOT", "f0_semitones"),
        category = "category",
        response = "response_category",
        group = "group",
        control = control_staninput(transform_type = "identity"),
        stanmodel = "NIX_ideal_adaptor"
      ),
      "requires exactly one cue"
    )
  }

  for (n_obs_exposure in c(0L, 1L, 3L)) {
    exp_1cue <- if (n_obs_exposure == 0L) {
      base_1cue$exposure[0L, ]
    } else {
      base_1cue$exposure[seq_len(n_obs_exposure), ]
    }
    exp_2cue <- if (n_obs_exposure == 0L) {
      base_2cue$exposure[0L, ]
    } else {
      base_2cue$exposure[seq_len(n_obs_exposure), ]
    }

    mnix_input <- new_ideal_adaptor_stanfit_input(
      exposure = exp_2cue,
      test = test_2cue,
      cues = c("VOT", "f0_semitones"),
      category = "category",
      response = "response_category",
      group = "group",
      control = control_staninput(transform_type = "identity"),
      stanmodel = "MNIX_ideal_adaptor"
    )
    expect_staninput_structure(
      mnix_input,
      expected_class = MNIX_IdealAdaptorStaninput,
      required_names = c(
        "N_exposure", "N_test", "x_mean_exposure", "x_ss_exposure", "p_cat"
      ),
      forbidden_names = c("x_sd_exposure")
    )

    expect_error(
      new_ideal_adaptor_stanfit_input(
        exposure = exp_1cue,
        test = test_1cue,
        cues = "VOT",
        category = "category",
        response = "response_category",
        group = "group",
        control = control_staninput(transform_type = "identity"),
        stanmodel = "MNIX_ideal_adaptor"
      ),
      "requires at least two cues"
    )
  }

  cue_combos <- list(
    list(cues = "VOT", base = base_1cue, test = test_1cue),
    list(cues = c("VOT", "f0_semitones"), base = base_2cue, test = test_2cue),
    list(
      cues = c("VOT", "f0_semitones", "vowel_duration"),
      base = base_3cue,
      test = test_3cue
    )
  )

  for (n_obs_exposure in c(0L, 1L, 3L)) {
    for (combo in cue_combos) {
      exp_data <- if (n_obs_exposure == 0L) {
        combo$base$exposure[0L, ]
      } else {
        combo$base$exposure[seq_len(n_obs_exposure), ]
      }

      niw_input <- new_ideal_adaptor_stanfit_input(
        exposure = exp_data,
        test = combo$test,
        cues = combo$cues,
        category = "category",
        response = "response_category",
        group = "group",
        control = control_staninput(transform_type = "identity"),
        stanmodel = "NIW_ideal_adaptor"
      )
      expect_staninput_structure(
        niw_input,
        expected_class = NIW_IdealAdaptorStaninput,
        required_names = c(
          "N_exposure", "N_test", "x_mean_exposure", "x_ss_exposure"
        ),
        forbidden_names = c("x_sd_exposure")
      )
    }
  }
})

test_that(
  "NIX supports one cue, rejects 2+ cues, matches NIX Stan input",
  {
  base_1cue <- example_exposure_test_data(model_family = "NIX", n_cues = 1L)
  base_2cue <- example_exposure_test_data(model_family = "MNIX", n_cues = 2L)
  base_1cue$exposure$group <- "g1"
  base_1cue$test$group <- "g1"
  base_2cue$exposure$group <- "g1"
  base_2cue$test$group <- "g1"
  test_1cue <- base_1cue$test[c(1, 41), ]

  for (n_obs_exposure in c(0L, 1L, 3L)) {
    exp_data <- if (n_obs_exposure == 0L) {
      base_1cue$exposure[0L, ]
    } else {
      base_1cue$exposure[seq_len(n_obs_exposure), ]
    }

    input <- new_ideal_adaptor_stanfit_input(
      exposure = exp_data,
      test = test_1cue,
      cues = "VOT",
      category = "category",
      response = "response_category",
      group = "group",
      control = control_staninput(transform_type = "identity"),
      stanmodel = "NIX_ideal_adaptor"
    )

    expected_counts <- as.vector(
      table(factor(exp_data$category, levels = c("/d/", "/t/")))
    )

    expect_staninput_structure(
      input,
      expected_class = NIX_IdealAdaptorStaninput,
      required_names = c(
        "N_exposure", "N_test", "x_mean_exposure", "x_ss_exposure"
      ),
      forbidden_names = c("x_sd_exposure")
    )
    expect_equal(as.vector(input@staninput@values$N_exposure), expected_counts)
    expect_equal(input@staninput@values$N_test, 2L)
  }

  expect_error(
    new_ideal_adaptor_stanfit_input(
      exposure = base_2cue$exposure[1L, ],
      test = base_2cue$test[c(1, 41), ],
      cues = c("VOT", "f0_semitones"),
      category = "category",
      response = "response_category",
      group = "group",
      control = control_staninput(transform_type = "identity"),
      stanmodel = "NIX_ideal_adaptor"
    ),
    "requires exactly one cue"
  )
})

test_that(
  "MNIX supports 2+ cues, rejects 1 cue, matches MNIX Stan input",
  {
  base_1cue <- example_exposure_test_data(model_family = "NIX", n_cues = 1L)
  base_2cue <- example_exposure_test_data(model_family = "MNIX", n_cues = 2L)
  base_1cue$exposure$group <- "g1"
  base_1cue$test$group <- "g1"
  base_2cue$exposure$group <- "g1"
  base_2cue$test$group <- "g1"
  test_2cue <- base_2cue$test[c(1, 41), ]

  for (n_obs_exposure in c(0L, 1L, 3L)) {
    exp_data <- if (n_obs_exposure == 0L) {
      base_2cue$exposure[0L, ]
    } else {
      base_2cue$exposure[seq_len(n_obs_exposure), ]
    }

    input <- new_ideal_adaptor_stanfit_input(
      exposure = exp_data,
      test = test_2cue,
      cues = c("VOT", "f0_semitones"),
      category = "category",
      response = "response_category",
      group = "group",
      control = control_staninput(transform_type = "identity"),
      stanmodel = "MNIX_ideal_adaptor"
    )

    expected_counts <- as.vector(
      table(factor(exp_data$category, levels = c("/d/", "/t/")))
    )

    expect_staninput_structure(
      input,
      expected_class = MNIX_IdealAdaptorStaninput,
      required_names = c(
        "N_exposure", "N_test", "x_mean_exposure", "x_ss_exposure", "p_cat"
      ),
      forbidden_names = c("x_sd_exposure")
    )
    expect_equal(as.vector(input@staninput@values$N_exposure), expected_counts)
    expect_equal(input@staninput@values$N_test, 2L)
  }

  expect_error(
    new_ideal_adaptor_stanfit_input(
      exposure = base_1cue$exposure[1L, ],
      test = base_1cue$test[c(1, 41), ],
      cues = "VOT",
      category = "category",
      response = "response_category",
      group = "group",
      control = control_staninput(transform_type = "identity"),
      stanmodel = "MNIX_ideal_adaptor"
    ),
    "requires at least two cues"
  )
})

test_that(
  "NIW supports 1, 2, and 2+ cues and matches NIW Stan input structure",
  {
  base_1cue <- example_exposure_test_data(model_family = "NIX", n_cues = 1L)
  base_2cue <- example_exposure_test_data(model_family = "MNIX", n_cues = 2L)
  base_3cue <- example_exposure_test_data(model_family = "NIW", n_cues = 3L)
  base_1cue$exposure$group <- "g1"
  base_1cue$test$group <- "g1"
  base_2cue$exposure$group <- "g1"
  base_2cue$test$group <- "g1"
  base_3cue$exposure$group <- "g1"
  base_3cue$test$group <- "g1"

  cue_combos <- list(
    list(
      cues = "VOT",
      base = base_1cue,
      test = base_1cue$test[c(1, 41), ]
    ),
    list(
      cues = c("VOT", "f0_semitones"),
      base = base_2cue,
      test = base_2cue$test[c(1, 41), ]
    ),
    list(
      cues = c("VOT", "f0_semitones", "vowel_duration"),
      base = base_3cue,
      test = base_3cue$test[c(1, 41), ]
    )
  )

  for (combo in cue_combos) {
    for (n_obs_exposure in c(0L, 1L, 3L)) {
      exp_data <- if (n_obs_exposure == 0L) {
        combo$base$exposure[0L, ]
      } else {
        combo$base$exposure[seq_len(n_obs_exposure), ]
      }

      input <- new_ideal_adaptor_stanfit_input(
        exposure = exp_data,
        test = combo$test,
        cues = combo$cues,
        category = "category",
        response = "response_category",
        group = "group",
        control = control_staninput(transform_type = "identity"),
        stanmodel = "NIW_ideal_adaptor"
      )

      expected_counts <- as.vector(
        table(factor(exp_data$category, levels = c("/d/", "/t/")))
      )

      expect_staninput_structure(
        input,
        expected_class = NIW_IdealAdaptorStaninput,
        required_names = c(
          "N_exposure", "N_test", "x_mean_exposure", "x_ss_exposure"
        ),
        forbidden_names = c("x_sd_exposure")
      )
      expect_equal(
        as.vector(input@staninput@values$N_exposure),
        expected_counts
      )
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
    response_category = factor(c("A", "B", "A")),
    cue1 = c(1.1, 2.1, 5.1),
    cue2 = c(3.1, 4.1, 7.1)
  )

  input1 <- new_ideal_adaptor_stanfit_input(
    exposure = exp1,
    test = test1,
    cues = c("cue1", "cue2"),
    category = "category",
    response = "response_category",
    group = "group",
    group_unique = "Condition",
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
      response = "response_category",
      group = "group",
      group_unique = "Condition",
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
    response = "response_category",
    group = "group",
    group_unique = "Condition",
    check_unique_group_identity = FALSE,
    control = control_staninput(transform_type = "identity"),
    stanmodel = "NIW_ideal_adaptor"
  )
  expect_equal(input_no_check@staninput@values$L, 2L)
  expect_equal(get_labels(input_no_check)$group, c("c1", "c2"))
})

test_that("IdealAdaptorStanfitInput stores standardized data and provides original variable names", {
  exp <- data.frame(
    MyCondition = factor("condA"),
    Subject = factor("s1"),
    TrueVowel = factor("vowel1"),
    f1 = 500,
    f2 = 1500
  )
  tst <- data.frame(
    MyCondition = factor("condA"),
    Subject = factor("s1"),
    UserChoice = factor("vowel1"),
    f1 = 520,
    f2 = 1480
  )

  input <- new_ideal_adaptor_stanfit_input(
    exposure = exp,
    test = tst,
    cues = c("f1", "f2"),
    category = "TrueVowel",
    response = "UserChoice",
    group = "Subject",
    group_unique = "MyCondition",
    control = control_staninput(transform_type = "identity"),
    stanmodel = "NIW_ideal_adaptor"
  )

  # Standard columns in @data
  expect_true(all(c("group", "group_unique", "category", "response_category", "f1", "f2", "Phase") %in% names(input@data)))
  expect_false("TrueVowel" %in% names(input@data))
  expect_false("UserChoice" %in% names(input@data))
  expect_false("Subject" %in% names(input@data))

  # get_original_variable_names
  orig_names <- get_original_variable_names(input)
  expect_equal(orig_names$group, "Subject")
  expect_equal(orig_names$group_unique, "MyCondition")
  expect_equal(orig_names$category, "TrueVowel")
  expect_equal(orig_names$response_category, "UserChoice")
  expect_equal(orig_names$cues, c("f1", "f2"))

  expect_equal(get_original_variable_names(input, "category"), "TrueVowel")
  expect_equal(get_original_variable_names(input, "response_category"), "UserChoice")
  expect_equal(get_original_variable_names(input, "group"), "Subject")
  expect_equal(get_original_variable_names(input, "group_unique"), "MyCondition")

  # get_data(original_names = FALSE)
  df_std <- get_data(input, original_names = FALSE)
  expect_true(all(c("group", "group_unique", "category", "response_category", "f1", "f2", "Phase") %in% names(df_std)))

  # get_data(original_names = TRUE)
  df_orig <- get_data(input, original_names = TRUE)
  expect_true(all(c("Subject", "MyCondition", "TrueVowel", "UserChoice", "f1", "f2", "Phase") %in% names(df_orig)))

  # get_exposure_data
  exp_std <- get_exposure_data(input, original_names = FALSE)
  expect_true("category" %in% names(exp_std))
  exp_orig <- get_exposure_data(input, original_names = TRUE)
  expect_true("TrueVowel" %in% names(exp_orig))

  # get_test_data
  tst_std <- get_test_data(input, original_names = FALSE)
  expect_true("response_category" %in% names(tst_std))
  tst_orig <- get_test_data(input, original_names = TRUE)
  expect_true("UserChoice" %in% names(tst_orig))
})

test_that("new_ideal_adaptor_stanfit_input handles group_unique and group.unique deprecation", {
  exp_df <- data.frame(
    Subject = c("s1", "s1", "s2", "s2", "s3", "s3", "s4", "s4"),
    MyCondition = c("cond1", "cond1", "cond1", "cond1", "cond2", "cond2", "cond2", "cond2"),
    TrueVowel = c("b", "d", "b", "d", "b", "d", "b", "d"),
    f1 = c(1, 2, 1, 2, 1, 2, 1, 2),
    f2 = c(3, 4, 3, 4, 3, 4, 3, 4),
    stringsAsFactors = FALSE
  )
  test_df <- data.frame(
    Subject = c("s1", "s2", "s3", "s4"),
    MyCondition = c("cond1", "cond1", "cond2", "cond2"),
    UserChoice = c("b", "d", "d", "b"),
    f1 = c(1.5, 2.5, 1.5, 2.5),
    f2 = c(3.5, 4.5, 3.5, 4.5),
    stringsAsFactors = FALSE
  )

  # Using group_unique
  input1 <- new_ideal_adaptor_stanfit_input(
    exposure = exp_df,
    test = test_df,
    cues = c("f1", "f2"),
    category = "TrueVowel",
    response = "UserChoice",
    group = "Subject",
    group_unique = "MyCondition"
  )
  expect_true(S7::S7_inherits(input1, IdealAdaptorStanfitInput))
  expect_equal(get_original_variable_names(input1, "group_unique"), "MyCondition")

  # Using deprecated group.unique emits warning
  expect_warning(
    input2 <- new_ideal_adaptor_stanfit_input(
      exposure = exp_df,
      test = test_df,
      cues = c("f1", "f2"),
      category = "TrueVowel",
      response = "UserChoice",
      group = "Subject",
      group.unique = "MyCondition"
    ),
    "deprecated"
  )
  expect_equal(get_original_variable_names(input2, "group_unique"), "MyCondition")
})

test_that("get_data, get_exposure_data, and get_test_data filtering and subsampling work", {
  exp_df <- data.frame(
    Subject = c("s1", "s1", "s2", "s2", "s3", "s3"),
    MyCondition = c("c1", "c1", "c1", "c1", "c2", "c2"),
    TrueVowel = c("b", "d", "b", "d", "b", "d"),
    f1 = c(1, 2, 1, 2, 1, 2),
    f2 = c(3, 4, 3, 4, 3, 4),
    stringsAsFactors = FALSE
  )
  test_df <- data.frame(
    Subject = c("s1", "s1", "s2", "s2", "s3", "s3"),
    MyCondition = c("c1", "c1", "c1", "c1", "c2", "c2"),
    UserChoice = c("b", "d", "b", "d", "d", "b"),
    f1 = c(1.5, 2.5, 1.5, 2.5, 1.5, 2.5),
    f2 = c(3.5, 4.5, 3.5, 4.5, 3.5, 4.5),
    stringsAsFactors = FALSE
  )

  input <- new_ideal_adaptor_stanfit_input(
    exposure = exp_df,
    test = test_df,
    cues = c("f1", "f2"),
    category = "TrueVowel",
    response = "UserChoice",
    group = "Subject",
    group_unique = "MyCondition"
  )

  # get_exposure_data: filter by groups
  exp_c1 <- get_exposure_data(input, groups = "c1")
  expect_equal(nrow(exp_c1), 4L)
  expect_true(all(exp_c1$group_unique == "c1"))

  # get_exposure_data: filter by invalid group errors
  expect_error(get_exposure_data(input, groups = "c99"), "Requested group.*not found")

  # get_exposure_data: filter by categories
  exp_b <- get_exposure_data(input, categories = "b")
  expect_equal(nrow(exp_b), 3L)
  expect_true(all(exp_b$category == "b"))

  # get_exposure_data: filter by invalid category errors
  expect_error(get_exposure_data(input, categories = "invalid_cat"), "Requested category.*not found")

  # get_exposure_data: subsample with n_samples
  exp_sub <- get_exposure_data(input, n_samples = 2)
  expect_equal(nrow(exp_sub), 2L)
  expect_error(get_exposure_data(input, n_samples = 0), "positive integer")

  # get_test_data: filter by groups and response_categories
  test_c2_d <- get_test_data(input, groups = "c2", response_categories = "d")
  expect_equal(nrow(test_c2_d), 1L)
  expect_equal(test_c2_d$response_category, "d")

  # get_test_data: filter by invalid response_categories errors
  expect_error(get_test_data(input, response_categories = "xyz"), "Requested response_category.*not found")

  # get_test_data: subsample with n_samples
  test_sub <- get_test_data(input, n_samples = 3)
  expect_equal(nrow(test_sub), 3L)

  # get_data: filter across phases
  all_c1 <- get_data(input, groups = "c1")
  expect_equal(nrow(all_c1), 8L) # 4 exposure + 4 test
  expect_true(all(all_c1$group_unique == "c1"))

  # get_data with categories (filters exposure) and response_categories
  # (filters test)
  all_filtered <- get_data(input, categories = "b", response_categories = "d")
  expect_true(
    all(all_filtered$category[all_filtered$Phase == "exposure"] == "b")
  )
  expect_true(
    all(all_filtered$response_category[all_filtered$Phase == "test"] == "d")
  )

  # get_data subsample
  all_sub <- get_data(input, n_samples = 4)
  expect_equal(nrow(all_sub), 4L)
})

test_that("get_data filters error on non-existing entities but message on zero-data", {
  exp_df <- data.frame(
    group = c("g1", "g1"),
    category = c("b", "d"),
    cue1 = c(1, 2),
    stringsAsFactors = FALSE
  )
  test_df <- data.frame(
    group = c("g1", "g1", "g_unexposed", "g_unexposed"),
    response_category = c("b", "d", "b", "b"),
    cue1 = c(1.5, 2.5, 1.2, 1.8),
    stringsAsFactors = FALSE
  )

  input <- new_ideal_adaptor_stanfit_input(
    exposure = exp_df,
    test = test_df,
    cues = "cue1",
    category = "category",
    response = "response_category",
    group = "group"
  )

  # 1. Non-existing group throws error
  expect_error(
    get_exposure_data(input, groups = "non_existent"),
    "Requested group.*not found"
  )
  expect_error(
    get_test_data(input, groups = "non_existent"),
    "Requested group.*not found"
  )
  expect_error(
    summary(input, groups = "non_existent"),
    "Requested group.*not found"
  )

  # 2. Existing group with zero exposure data emits message, not error
  expect_message(
    res_zero <- get_exposure_data(input, groups = "g_unexposed"),
    "No exposure data found for group"
  )
  expect_equal(nrow(res_zero), 0L)

  # 3. Requesting both groups emits message for g_unexposed and returns g1 data
  expect_message(
    res_both <- get_exposure_data(input, groups = c("g1", "g_unexposed")),
    "No exposure data found for group"
  )
  expect_equal(nrow(res_both), 2L)
  expect_true(all(res_both$group == "g1"))

  # 4. Non-existing category throws error
  expect_error(
    get_exposure_data(input, categories = "non_existent_cat"),
    "Requested category.*not found"
  )
  expect_error(
    summary(input, categories = "non_existent_cat"),
    "Requested category.*not found"
  )

  # 5. Non-existing response_category throws error
  expect_error(
    get_test_data(input, response_categories = "non_existent_resp"),
    "Requested response_category.*not found"
  )

  # 6. Existing response_category with zero data in filtered group messages
  expect_message(
    res_rcat <- get_test_data(
      input,
      groups = "g_unexposed",
      response_categories = "d"
    ),
    "No test data found for response_category"
  )
  expect_equal(nrow(res_rcat), 0L)

  # 7. summary on input with unexposed group summarizes existing data without error
  summ <- summary(input)
  expect_true(S7::S7_inherits(summ, Summary_IdealAdaptorStanfitInput))
  expect_true("g1" %in% summ@exposure_statistics$group)
  expect_false("g_unexposed" %in% summ@exposure_statistics$group)
  expect_equal(summ@test_summary$n_observations, 4L)

  # 8. summary with explicit group filter
  summ_g1 <- summary(input, groups = "g1")
  expect_equal(summ_g1@groups, "g1")
  expect_equal(summ_g1@test_summary$n_observations, 2L)
})

test_that(
  "constructor-generated data are compatible with matching Stan programs",
  {
    base_1cue <- example_exposure_test_data(model_family = "NIX", n_cues = 1L)
    base_2cue <- example_exposure_test_data(model_family = "MNIX", n_cues = 2L)
    base_3cue <- example_exposure_test_data(model_family = "NIW", n_cues = 3L)
    base_1cue$exposure$group <- "g1"
    base_1cue$test$group <- "g1"
    base_2cue$exposure$group <- "g1"
    base_2cue$test$group <- "g1"
    base_3cue$exposure$group <- "g1"
    base_3cue$test$group <- "g1"

    test_1cue <- base_1cue$test[c(1, 41), ]
    test_2cue <- base_2cue$test[c(1, 41), ]
    test_3cue <- base_3cue$test[c(1, 41), ]

    for (n_obs_exposure in c(0L, 1L, 3L)) {
      exp_1cue <- if (n_obs_exposure == 0L) {
        base_1cue$exposure[0L, ]
      } else {
        base_1cue$exposure[seq_len(n_obs_exposure), ]
      }

      input <- new_ideal_adaptor_stanfit_input(
        exposure = exp_1cue,
        test = test_1cue,
        cues = "VOT",
        category = "category",
        response = "response_category",
        group = "group",
        control = control_staninput(transform_type = "identity"),
        stanmodel = "NIX_ideal_adaptor"
      )

      run_fixed_param_stan_program(
        stan_file = system.file(
          "stan", "NIX_ideal_adaptor.stan",
          package = "MVBeliefUpdatr"
        ),
        input = input,
        model_name = "NIX_compat"
      )
    }

    mnix_combos <- list(
      list(cues = c("VOT", "f0_semitones"), base = base_2cue, test = test_2cue),
      list(
        cues = c("VOT", "f0_semitones", "vowel_duration"),
        base = base_3cue,
        test = test_3cue
      )
    )

    for (n_obs_exposure in c(0L, 1L, 3L)) {
      for (combo in mnix_combos) {
        exp_mnix <- if (n_obs_exposure == 0L) {
          combo$base$exposure[0L, ]
        } else {
          combo$base$exposure[seq_len(n_obs_exposure), ]
        }

        input <- new_ideal_adaptor_stanfit_input(
          exposure = exp_mnix,
          test = combo$test,
          cues = combo$cues,
          category = "category",
          response = "response_category",
          group = "group",
          control = control_staninput(transform_type = "identity"),
          stanmodel = "MNIX_ideal_adaptor"
        )

        run_fixed_param_stan_program(
          stan_file = system.file(
            "stan", "MNIX_ideal_adaptor.stan",
            package = "MVBeliefUpdatr"
          ),
          input = input,
          model_name = "MNIX_compat"
        )
      }
    }

    niw_combos <- list(
      list(cues = "VOT", base = base_1cue, test = test_1cue),
      list(cues = c("VOT", "f0_semitones"), base = base_2cue, test = test_2cue),
      list(
        cues = c("VOT", "f0_semitones", "vowel_duration"),
        base = base_3cue,
        test = test_3cue
      )
    )

    for (n_obs_exposure in c(0L, 1L, 3L)) {
      for (combo in niw_combos) {
        exp_data <- if (n_obs_exposure == 0L) {
          combo$base$exposure[0L, ]
        } else {
          combo$base$exposure[seq_len(n_obs_exposure), ]
        }

        input <- new_ideal_adaptor_stanfit_input(
          exposure = exp_data,
          test = combo$test,
          cues = combo$cues,
          category = "category",
          response = "response_category",
          group = "group",
          control = control_staninput(transform_type = "identity"),
          stanmodel = "NIW_ideal_adaptor"
        )

        run_fixed_param_stan_program(
          stan_file = system.file(
            "stan", "NIW_ideal_adaptor.stan",
            package = "MVBeliefUpdatr"
          ),
          input = input,
          model_name = "NIW_compat"
        )
      }
    }
  }
)


