# test-38-S7-prior.R
# Tests for MVBU_Prior, prior(), and dual-space prior resolution

test_that("prior() constructs valid MVBU_Prior objects", {
  p_tau <- prior("cauchy(0, 5)", class = "tau")
  expect_true(S7::S7_inherits(p_tau, MVBU_Prior))
  expect_equal(p_tau@declaration, "cauchy(0, 5)")
  expect_equal(p_tau@par_class, "tau")
  expect_equal(p_tau@par, "")
  expect_equal(p_tau@cue_space, "transformed")

  # Numeric shorthand for tau
  p_num_tau <- prior(5, class = "tau")
  expect_equal(p_num_tau@declaration, "cauchy(0, 5)")

  # Numeric shorthand for cor
  p_num_cor <- prior(2, class = "cor")
  expect_equal(p_num_cor@declaration, "lkj(2)")

  # Original cue space specification
  p_orig <- prior("cauchy(0, 100)", class = "tau", par = "f1",
                  cue_space = "original")
  expect_equal(p_orig@cue_space, "original")
  expect_equal(p_orig@par, "f1")
})

test_that("prior() input validation catches invalid declarations and classes", {
  expect_error(prior(), "Argument `declaration` must be provided")
  expect_error(prior("", class = "tau"), "non-empty character string")
  expect_error(prior("cauchy(0, 5)", class = "unknown_class"), "arg")
  expect_error(
    prior("cauchy(0, 5)", class = "tau", cue_space = "invalid_space"),
    "arg"
  )
})

test_that("MVBU_Prior S3 methods work as expected", {
  p1 <- prior("cauchy(0, 5)", class = "tau")
  p2 <- prior("lkj(1)", class = "cor")
  p3 <- prior("normal(0, 1)", class = "m")

  # Combine via c()
  combined <- c(p1, p2, p3)
  expect_true(S7::S7_inherits(combined, MVBU_Prior))
  expect_equal(length(combined), 3L)

  # Nested list combine
  combined_list <- c(p1, list(p2, p3))
  expect_equal(length(combined_list), 3L)

  # Subsetting
  sub_p <- combined[c(1, 3)]
  expect_equal(length(sub_p), 2L)
  expect_equal(sub_p@par_class, c("tau", "m"))

  # as.data.frame
  df <- as.data.frame(combined)
  expect_s3_class(df, "data.frame")
  expect_equal(
    colnames(df),
    c("declaration", "class", "category", "cue", "par", "cue_space")
  )
  expect_equal(df$class, c("tau", "cor", "m"))

  # as_tibble
  tbl <- tibble::as_tibble(combined)
  expect_s3_class(tbl, "tbl_df")
  expect_equal(nrow(tbl), 3L)

  # print output
  out <- capture.output(print(combined))
  expect_true(any(grepl("cauchy\\(0, 5\\)", out)))
  expect_true(any(grepl("lkj\\(1\\)", out)))
})

test_that(".extract_tau_scale and .extract_lkj_eta parse values correctly", {
  expect_equal(MVBeliefUpdatr:::.extract_tau_scale("cauchy(0, 5)"), 5)
  expect_equal(MVBeliefUpdatr:::.extract_tau_scale("cauchy(0, 2.5)"), 2.5)
  expect_equal(MVBeliefUpdatr:::.extract_tau_scale("5"), 5)
  expect_error(
    MVBeliefUpdatr:::.extract_tau_scale("invalid_declaration"),
    "Could not extract a positive scale parameter"
  )

  expect_equal(MVBeliefUpdatr:::.extract_lkj_eta("lkj(1)"), 1)
  expect_equal(MVBeliefUpdatr:::.extract_lkj_eta("lkj_corr_cholesky(2.5)"), 2.5)
  expect_equal(MVBeliefUpdatr:::.extract_lkj_eta("1"), 1)
  expect_error(
    MVBeliefUpdatr:::.extract_lkj_eta("invalid_declaration"),
    "Could not extract a positive eta parameter"
  )
})

test_that(".resolve_priors resolves default dual-space priors correctly", {
  # Setup mock affine transform
  # Cues: c0 (scale = 0.1, orig multiplier = 10)
  #       c1 (scale = 0.05, mult = 20)
  scale_mat <- diag(c(0.1, 0.05))
  shift_vec <- c(-50, -100)
  transform <- MVBU_TransformInformation(
    cue.labels = c("c0", "c1"),
    type = "standardize",
    shift = shift_vec,
    SCALE = scale_mat
  )

  res <- MVBeliefUpdatr:::.resolve_priors(
    priors = list(),
    stanmodel = "NIW_ideal_adaptor",
    cues = c("c0", "c1"),
    category_labels = c("cat1", "cat2"),
    transform = transform
  )

  expect_equal(res$tau_scale, c(5, 5))
  expect_equal(res$L_omega_eta, 1)

  # Transformed space check
  trans_tbl <- as.data.frame(res$transformed)
  expect_true("cauchy(0, 5)" %in% trans_tbl$declaration)
  expect_true("lkj(1)" %in% trans_tbl$declaration)
  expect_true("normal(0, tau)" %in% trans_tbl$declaration)
  expect_true("dirichlet(1)" %in% trans_tbl$declaration)

  # Original space check
  orig_tbl <- as.data.frame(res$original)
  # c0 tau: 5 * 10 = 50; c1 tau: 5 * 20 = 100
  expect_true(
    any(orig_tbl$declaration == "cauchy(0, 50)" & orig_tbl$par == "c0")
  )
  expect_true(
    any(orig_tbl$declaration == "cauchy(0, 100)" & orig_tbl$par == "c1")
  )
  expect_true(any(orig_tbl$declaration == "lkj(1)" & orig_tbl$class == "cor"))
  # m in original space: shift is -mean, so center is 50 for c0 and 100 for c1
  expect_true(
    any(orig_tbl$declaration == "normal(50, tau)" & orig_tbl$par == "c0")
  )
  expect_true(
    any(orig_tbl$declaration == "normal(100, tau)" & orig_tbl$par == "c1")
  )
})

test_that(
  ".resolve_priors supports custom priors in transformed and original space",
  {
  scale_mat <- diag(c(0.1, 0.05))
  shift_vec <- c(-50, -100)
  transform <- MVBU_TransformInformation(
    cue.labels = c("c0", "c1"),
    type = "standardize",
    shift = shift_vec,
    SCALE = scale_mat
  )

  # User specifies transformed tau and LKJ
  priors_trans <- c(
    prior("cauchy(0, 2.5)", class = "tau"),
    prior("lkj(2)", class = "cor")
  )
  res_trans <- MVBeliefUpdatr:::.resolve_priors(
    priors = priors_trans,
    stanmodel = "NIW_ideal_adaptor",
    cues = c("c0", "c1"),
    category_labels = c("cat1", "cat2"),
    transform = transform
  )
  expect_equal(res_trans$tau_scale, c(2.5, 2.5))
  expect_equal(res_trans$L_omega_eta, 2)
  orig_df <- as.data.frame(res_trans$original)
  expect_true(any(orig_df$declaration == "cauchy(0, 25)" & orig_df$par == "c0"))
  expect_true(any(orig_df$declaration == "cauchy(0, 50)" & orig_df$par == "c1"))

  # User specifies cue-specific original space prior for c0
  priors_orig <- prior("cauchy(0, 100)", class = "tau", par = "c0",
                       cue_space = "original")
  res_orig <- MVBeliefUpdatr:::.resolve_priors(
    priors = priors_orig,
    stanmodel = "NIW_ideal_adaptor",
    cues = c("c0", "c1"),
    category_labels = c("cat1", "cat2"),
    transform = transform
  )
  # c0: 100 original / 10 mult = 10 transformed. c1: default 5
  expect_equal(res_orig$tau_scale, c(10, 5))
})

test_that("control_staninput deprecation warning triggers appropriately", {
  # Default call triggers NO warning
  expect_no_warning(ctrl <- control_staninput())
  expect_null(ctrl$tau_scale)
  expect_null(ctrl$L_omega_eta)

  # Specifying tau_scale triggers deprecation warning
  expect_warning(
    ctrl_dep <- control_staninput(tau_scale = 10),
    "deprecated"
  )
  expect_equal(ctrl_dep$tau_scale, 10)

  # Specifying L_omega_eta triggers deprecation warning
  expect_warning(
    ctrl_dep2 <- control_staninput(L_omega_eta = 2),
    "deprecated"
  )
  expect_equal(ctrl_dep2$L_omega_eta, 2)
})

test_that("get_stan_priors works for MVBU_Prior objects", {
  p_trans <- prior("cauchy(0, 5)", class = "tau")
  p_orig <- prior("cauchy(0, 5)", class = "tau", cue_space = "original")
  p_both <- c(p_trans, p_orig)
  expect_identical(get_stan_priors(p_trans, space = "transformed"), p_trans)
  expect_identical(get_stan_priors(p_trans, space = "original"), MVBU_Prior())
  expect_identical(get_stan_priors(p_both, space = "transformed"), p_trans)
  expect_identical(get_stan_priors(p_both, space = "original"), p_orig)
})
