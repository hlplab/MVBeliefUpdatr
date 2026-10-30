# test-39-S7-fixed.R
# Tests for MVBU_FixedParameter, fix_parameter(), and fixed parameter resolution

test_that("fix_parameter() constructs valid MVBU_FixedParameter objects", {
  # m_0
  f_m <- fix_parameter(c(0, 0), class = "m_0", category = "b")
  expect_true(S7::S7_inherits(f_m, MVBU_FixedParameter))
  expect_equal(f_m@par_class, "m_0")
  expect_equal(f_m@category, "b")
  expect_equal(f_m@cue_space, "original")
  expect_equal(f_m@value[[1]], c(0, 0))

  # expected_mu_0 canonicalizes to m_0
  f_emu <- fix_parameter(c(1, 2), class = "expected_mu_0", category = "d")
  expect_equal(f_emu@par_class, "m_0")
  expect_equal(f_emu@category, "d")

  # S_0 and expected_Sigma_0
  Sigma_mat <- diag(2)
  f_s <- fix_parameter(Sigma_mat, class = "S_0", category = "b")
  expect_equal(f_s@par_class, "S_0")
  expect_equal(f_s@category, "b")

  f_es <- fix_parameter(Sigma_mat, class = "expected_Sigma_0", category = "b")
  expect_equal(f_es@par_class, "expected_Sigma_0")
  expect_equal(f_es@category, "b")

  # kappa_0 and nu_0
  f_k <- fix_parameter(5, class = "kappa_0", category = "b")
  expect_equal(f_k@par_class, "kappa_0")
  f_nu <- fix_parameter(10, class = "nu_0")
  expect_equal(f_nu@par_class, "nu_0")
  expect_equal(f_nu@category, "")

  # Global parameters: lapse_rate, lapse_bias, category_prior
  f_lr <- fix_parameter(0.02, class = "lapse_rate")
  expect_equal(f_lr@par_class, "lapse_rate")
  expect_equal(f_lr@category, "")

  f_lb <- fix_parameter(c(0.5, 0.5), class = "lapse_bias")
  expect_equal(f_lb@par_class, "lapse_bias")

  f_cp <- fix_parameter(c(0.5, 0.5), class = "category_prior")
  expect_equal(f_cp@par_class, "category_prior")

  # mu_0 and Sigma_0 (Stage 6)
  f_mu0 <- fix_parameter(c(1, 2), class = "mu_0", category = "b")
  expect_equal(f_mu0@par_class, "mu_0")
  expect_equal(f_mu0@category, "b")
  expect_equal(f_mu0@cue_space, "original")
  expect_equal(f_mu0@value[[1]], c(1, 2))

  f_sig0 <- fix_parameter(Sigma_mat, class = "Sigma_0", category = "b")
  expect_equal(f_sig0@par_class, "Sigma_0")
  expect_equal(f_sig0@category, "b")

  # weights
  f_w <- fix_parameter(1.5, class = "weights", cue = "f1")
  expect_equal(f_w@par_class, "weights")
  expect_equal(f_w@cue, "f1")

  # String assignment syntax
  f_str <- fix_parameter("lapse_rate = 0.05")
  expect_equal(f_str@par_class, "lapse_rate")
  expect_equal(f_str@value[[1]], 0.05)
})

test_that(
  "fix_parameter() validates constraints and throws informative errors",
  {
    # Value required
    expect_error(fix_parameter(), "must be provided")

  # Global parameters cannot be category-specific
  expect_error(
    fix_parameter(0.05, class = "lapse_rate", category = "b"),
    "is a global parameter and cannot be category-specific"
  )
  expect_error(
    fix_parameter(c(0.5, 0.5), class = "category_prior", category = "b"),
    "is a global parameter and cannot be category-specific"
  )
  expect_error(
    fix_parameter(c(0.5, 0.5), class = "lapse_bias", category = "b"),
    "is a global parameter and cannot be category-specific"
  )

  # Cue specificity constraints
  expect_error(
    fix_parameter(0.05, class = "lapse_rate", cue = "f1"),
    "cannot be cue-specific"
  )
  expect_error(
    fix_parameter(5, class = "kappa_0", cue = "f1"),
    "cannot be cue-specific"
  )
  expect_error(
    fix_parameter(10, class = "nu_0", cue = "f1"),
    "cannot be cue-specific"
  )
})

test_that("MVBU_FixedParameter S3 methods work as expected", {
  f1 <- fix_parameter(0.05, class = "lapse_rate")
  f2 <- fix_parameter(c(0, 0), class = "m_0", category = "b")
  f3 <- fix_parameter(diag(2), class = "S_0", category = "b")

  # Combine via c()
  combined <- c(f1, f2, f3)
  expect_true(S7::S7_inherits(combined, MVBU_FixedParameter))
  expect_equal(length(combined), 3L)

  # names()
  expect_equal(names(combined), c("lapse_rate", "m_0", "S_0"))

  # Subsetting
  sub_f <- combined[c(1, 3)]
  expect_equal(length(sub_f), 2L)
  expect_equal(sub_f@par_class, c("lapse_rate", "S_0"))

  # as.data.frame
  df <- as.data.frame(combined)
  expect_s3_class(df, "data.frame")
  expect_equal(
    colnames(df),
    c("class", "category", "cue", "cue_space", "value")
  )
  expect_equal(df$class, c("lapse_rate", "m_0", "S_0"))

  # as_tibble
  tbl <- tibble::as_tibble(combined)
  expect_s3_class(tbl, "tbl_df")
  expect_equal(nrow(tbl), 3L)

  # print output
  out <- capture.output(print(combined))
  expect_true(any(grepl("lapse_rate", out)))
  expect_true(any(grepl("m_0", out)))
  expect_true(any(grepl("S_0", out)))
})

test_that(
  ".as_fixed_parameters() rejects legacy lists and accepts MVBU_FixedParameter",
  {
    # NULL returns empty MVBU_FixedParameter
    empty_fp <- MVBeliefUpdatr:::.as_fixed_parameters(NULL)
    expect_true(S7::S7_inherits(empty_fp, MVBU_FixedParameter))
    expect_equal(length(empty_fp), 0L)

    # MVBU_FixedParameter returned as-is
    fp <- fix_parameter(0.05, class = "lapse_rate")
    expect_identical(MVBeliefUpdatr:::.as_fixed_parameters(fp), fp)

    # Legacy list is rejected without silent coercion
    expect_error(
      MVBeliefUpdatr:::.as_fixed_parameters(list(lapse_rate = 0.05)),
      "must be an MVBU_FixedParameter object"
    )
  }
)

test_that(
  ".resolve_fixed_parameters() correctly resolves and transforms values",
  {
    transform <- MVBeliefUpdatr::MVBU_TransformInformation(
      cue.labels = c("f1", "f2"),
      type = "standardize",
      shift = c(10, 20),
      SCALE = matrix(c(2, 0, 0, 3), nrow = 2, byrow = TRUE)
    )
    categories <- c("b", "d")

    fp <- c(
      fix_parameter(0.03, class = "lapse_rate"),
      fix_parameter(c(0.6, 0.4), class = "category_prior"),
      fix_parameter(
        c(10, 20),
        class = "m_0",
        category = "b",
        cue_space = "original"
      ),
      fix_parameter(
        matrix(c(1, 0, 0, 1), 2, 2),
        class = "S_0",
        category = "b",
        cue_space = "original"
      ),
      fix_parameter(
        matrix(c(2, 0, 0, 2), 2, 2),
        class = "expected_Sigma_0",
        category = "d",
        cue_space = "original"
      ),
      fix_parameter(5, class = "kappa_0", category = "b"),
      fix_parameter(8, class = "nu_0", category = "b")
    )

    resolved <- MVBeliefUpdatr:::.resolve_fixed_parameters(
      fixed = fp,
      categories = categories,
      cues = c("f1", "f2"),
      transform = transform
    )

    expect_equal(resolved$lapse_rate$known, 1L)
    expect_equal(as.numeric(resolved$lapse_rate$data), 0.03)

    expect_equal(resolved$category_prior$known, 1L)
    expect_equal(as.numeric(resolved$category_prior$data), c(0.6, 0.4))

    # m_0 transformed: SCALE %*% (val - shift) = SCALE %*% (0) = c(0, 0)
    expect_equal(resolved$m_0$M_known, 1L)
    expect_equal(resolved$m_0$idx_known, as.array(1L))
    expect_equal(as.numeric(resolved$m_0$data[1, ]), c(0, 0))

    # S_0 (mode 1) transformed: SCALE %*% S %*% t(SCALE) = diag(4, 9)
    expect_equal(resolved$S_0$M_known, 2L)
    expect_equal(resolved$S_0$mode_known, as.array(c(1L, 2L)))
    expect_equal(resolved$S_0$data[1, , ], matrix(c(4, 0, 0, 9), 2, 2))
    # expected_Sigma_0 (mode 2) transformed: diag(4, 9) * 2 = diag(8, 18)
    expect_equal(resolved$S_0$data[2, , ], matrix(c(8, 0, 0, 18), 2, 2))

    # kappa_0 and nu_0
    expect_equal(resolved$kappa_0$M_known, 1L)
    expect_equal(as.numeric(resolved$kappa_0$data[1]), 5)
    expect_equal(resolved$nu_0$M_known, 1L)
    expect_equal(as.numeric(resolved$nu_0$data[1]), 8)
  }
)

test_that(
  ".resolve_fixed_parameters() supports known mu_0 and known Sigma_0",
  {
    categories <- c("b", "d")
    transform <- MVBeliefUpdatr::MVBU_TransformInformation(
      cue.labels = c("f1", "f2"),
      type = "standardize",
      shift = c(10, 20),
      SCALE = matrix(c(2, 0, 0, 3), nrow = 2, byrow = TRUE)
    )

    # 1. Valid known mu_0 and known Sigma_0
    fp <- list(
      fix_parameter(c(10, 20), class = "mu_0", category = "b"),
      fix_parameter(
        matrix(c(1, 0, 0, 1), 2, 2),
        class = "Sigma_0",
        category = "d"
      )
    )

    resolved <- MVBeliefUpdatr:::.resolve_fixed_parameters(
      fixed = fp,
      categories = categories,
      cues = c("f1", "f2"),
      transform = transform
    )

    expect_equal(resolved$mu_0$is_known, as.array(c(1L, 0L)))
    expect_equal(as.numeric(resolved$mu_0$data[1, ]), c(0, 0))
    expect_equal(resolved$Sigma_0$is_known, as.array(c(0L, 1L)))
    expect_equal(resolved$Sigma_0$data[2, , ], matrix(c(4, 0, 0, 9), 2, 2))

    # Known mu_0 is added to m_0$idx_known to prevent Stan from sampling it
    expect_true(1L %in% resolved$m_0$idx_known)
    # Known Sigma_0 is added to S_0$idx_known with mode 1
    expect_true(2L %in% resolved$S_0$idx_known)

    # 2. Conflicting specification: both m_0 and mu_0 for same category
    fp_conflict_m <- list(
      fix_parameter(c(0, 0), class = "m_0", category = "b"),
      fix_parameter(c(10, 20), class = "mu_0", category = "b")
    )
    expect_error(
      MVBeliefUpdatr:::.resolve_fixed_parameters(
        fixed = fp_conflict_m,
        categories = categories,
        cues = c("f1", "f2"),
        transform = transform
      ),
      "Cannot fix both 'm_0' and known 'mu_0' for category 'b'"
    )

    # 3. Conflicting specification: both S_0 and Sigma_0 for same category
    fp_conflict_s <- list(
      fix_parameter(diag(2), class = "S_0", category = "d"),
      fix_parameter(diag(2), class = "Sigma_0", category = "d")
    )
    expect_error(
      MVBeliefUpdatr:::.resolve_fixed_parameters(
        fixed = fp_conflict_s,
        categories = categories,
        cues = c("f1", "f2"),
        transform = transform
      ),
      "Cannot fix both 'S_0'/'expected_Sigma_0' and known 'Sigma_0'"
    )

    # 4. Non-positive-definite Sigma_0 errors
    fp_bad_sigma <- list(
      fix_parameter(
        matrix(c(1, 2, 2, 1), 2, 2),
        class = "Sigma_0",
        category = "b"
      )
    )
    expect_error(
      MVBeliefUpdatr:::.resolve_fixed_parameters(
        fixed = fp_bad_sigma,
        categories = categories,
        cues = c("f1", "f2"),
        transform = transform
      ),
      "must be positive-definite"
    )
  }
)

