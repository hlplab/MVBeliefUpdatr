test_that("update_template performs an S7 NIW batch update", {
  model <- example_niw_ideal_adaptor(n_cues = 2)
  categories <- get_category_labels(model)[1:2]
  exposure <- tibble::tibble(
    category = rep(categories, each = 2),
    VOT = c(10, 20, 40, 50),
    f0_semitones = c(-1, 0, 1, 2)
  )

  updated <- update_template(model, exposure)
  before <- model@category_template@representations
  after <- updated@category_template@representations

  expect_true(S7::S7_inherits(updated, NIW_IdealAdaptor))
  expect_equal(after[[1]]@kappa, before[[1]]@kappa + 2)
  expect_equal(after[[1]]@nu, before[[1]]@nu + 2)
  expect_false(identical(after[[1]]@m, before[[1]]@m))
  expect_false(identical(after[[1]]@S, before[[1]]@S))
  expect_equal(after[[2]]@kappa, before[[2]]@kappa + 2)
  expect_equal(after[[2]]@nu, before[[2]]@nu + 2)
})

test_that("update_category_representation updates an NIW representation", {
  representation <- example_niw_category_representation(n_cues = 2, nu = 30)
  updated <- update_category_representation(
    representation, x_N = 2, x_mean = c(15, .5),
    x_css = matrix(c(50, 5, 5, .5), nrow = 2)
  )
  expect_true(S7::S7_inherits(updated, NIW_CategoryRepresentation))
  expect_equal(updated@kappa, representation@kappa + 2)
  expect_equal(updated@nu, representation@nu + 2)
  expect_false(identical(updated@m, representation@m))
  expect_false(identical(updated@S, representation@S))
})

test_that("update_template preserves categories without exposure", {
  model <- example_niw_ideal_adaptor(n_cues = 1)
  category <- get_category_labels(model)[1]
  exposure <- tibble::tibble(category = category, VOT = c(10, 20))
  updated <- update_template(model, exposure)
  expect_identical(
    updated@category_template@representations[[2]]@m,
    model@category_template@representations[[2]]@m
  )
  expect_identical(
    updated@category_template@representations[[2]]@S,
    model@category_template@representations[[2]]@S
  )
})

test_that("update_template rejects incompatible S7 update data", {
  model <- example_niw_ideal_adaptor(n_cues = 1)
  expect_error(update_template(model, data.frame(category = "A", wrong_cue = 1)))
  expect_error(update_template(model, data.frame(wrong_category = "A", VOT = 1)))
})

test_that("update_template supports incremental history and update methods", {
  model <- example_niw_ideal_adaptor(n_cues = 1)
  observations <- tibble::tibble(
    category = get_category_labels(model)[1],
    VOT = 10
  )

  history <- update_template(
    model,
    observations,
    updating = "incremental",
    keep_history = TRUE,
    lapse_treatment = "no_lapses",
    noise_treatment = "no_noise",
    update_method = "label-certain"
  )
  expect_type(history, "list")
  expect_length(history, 2)
  expect_true(all(vapply(history, function(x) S7::S7_inherits(x, NIW_IdealAdaptor), logical(1))))

  unlabeled <- observations[0, c("VOT"), drop = FALSE]
  expect_true(S7::S7_inherits(
    update_template(model, unlabeled, updating = "incremental", update_method = "nolabel-uniform"),
    NIW_IdealAdaptor
  ))
  expect_error(update_template(model, observations, updating = "batch", update_method = "nolabel-uniform"))
  expect_error(update_template(model, observations, update_method = "not-a-method"))
})

test_that("update_template performs an S7 NIX batch and incremental update", {
  model <- example_nix_ideal_adaptor()
  categories <- get_category_labels(model)[1:2]
  exposure <- tibble::tibble(
    category = rep(categories, each = 2),
    VOT = c(10, 20, 40, 50)
  )

  # Batch update
  updated <- update_template(model, exposure, updating = "batch")
  before <- model@category_template@representations
  after <- updated@category_template@representations

  expect_true(S7::S7_inherits(updated, NIX_IdealAdaptor))
  expect_equal(after[[1]]@kappa, before[[1]]@kappa + 2)
  expect_equal(after[[1]]@nu, before[[1]]@nu + 2)
  expect_false(identical(after[[1]]@m, before[[1]]@m))
  expect_false(identical(after[[1]]@sigma2, before[[1]]@sigma2))
  expect_equal(after[[2]]@kappa, before[[2]]@kappa + 2)
  expect_equal(after[[2]]@nu, before[[2]]@nu + 2)

  # Incremental update with history
  history <- update_template(model, exposure, updating = "incremental", keep_history = TRUE)
  expect_type(history, "list")
  expect_length(history, nrow(exposure) + 1L)
  expect_true(all(vapply(history, function(m) S7::S7_inherits(m, NIX_IdealAdaptor), logical(1))))

  # Final incremental model matches batch update parameters
  final_inc <- history[[length(history)]]
  expect_equal(final_inc@category_template@representations[[1]]@m, after[[1]]@m)
  expect_equal(final_inc@category_template@representations[[1]]@sigma2, after[[1]]@sigma2)
  expect_equal(final_inc@category_template@representations[[1]]@kappa, after[[1]]@kappa)
  expect_equal(final_inc@category_template@representations[[1]]@nu, after[[1]]@nu)
})

test_that("update_category_representation updates an NIX representation", {
  rep0 <- example_nix_category_representation()
  updated <- update_category_representation(
    rep0, x_N = 3, x_mean = 25, x_css = 12
  )
  expect_true(S7::S7_inherits(updated, NIX_CategoryRepresentation))
  expect_equal(updated@kappa, rep0@kappa + 3)
  expect_equal(updated@nu, rep0@nu + 3)
  expect_false(identical(updated@m, rep0@m))
  expect_false(identical(updated@sigma2, rep0@sigma2))
})

test_that("update_template performs an S7 MNIX batch and incremental update", {
  model <- example_mnix_ideal_adaptor()
  categories <- get_category_labels(model)[1:2]
  exposure <- tibble::tibble(
    category = rep(categories, each = 2),
    VOT = c(10, 20, 40, 50),
    f0_semitones = c(-1, 0, 1, 2)
  )

  # Batch update
  updated <- update_template(model, exposure, updating = "batch")
  before <- model@category_template@representations
  after <- updated@category_template@representations

  expect_true(S7::S7_inherits(updated, MNIX_IdealAdaptor))
  expect_equal(after[[1]]@kappa, before[[1]]@kappa + 2)
  expect_equal(after[[1]]@nu, before[[1]]@nu + 2)
  expect_false(identical(after[[1]]@m, before[[1]]@m))
  expect_false(identical(after[[1]]@sigma2, before[[1]]@sigma2))

  # Incremental update with history
  history <- update_template(model, exposure, updating = "incremental", keep_history = TRUE)
  expect_type(history, "list")
  expect_length(history, nrow(exposure) + 1L)
  expect_true(all(vapply(history, function(m) S7::S7_inherits(m, MNIX_IdealAdaptor), logical(1))))

  # Final incremental model matches batch update parameters
  final_inc <- history[[length(history)]]
  expect_equal(final_inc@category_template@representations[[1]]@m, after[[1]]@m)
  expect_equal(final_inc@category_template@representations[[1]]@sigma2, after[[1]]@sigma2)
  expect_equal(final_inc@category_template@representations[[1]]@kappa, after[[1]]@kappa)
  expect_equal(final_inc@category_template@representations[[1]]@nu, after[[1]]@nu)
})

test_that("update_category_representation updates an MNIX representation", {
  rep0 <- example_mnix_category_representation()
  D <- length(get_cue_labels(rep0))
  updated <- update_category_representation(
    rep0, x_N = 4, x_mean = c(15, 0.5), x_css = diag(c(20, 1), D)
  )
  expect_true(S7::S7_inherits(updated, MNIX_CategoryRepresentation))
  expect_equal(updated@kappa, rep0@kappa + 4)
  expect_equal(updated@nu, rep0@nu + 4)
  expect_false(identical(updated@m, rep0@m))
  expect_false(identical(updated@sigma2, rep0@sigma2))
})

test_that("1-cue analytical parity between NIX and NIW forward updating", {
  # Setup matching 1-cue category representations:
  # In 1D, NIX has scale parameter sigma2, while NIW has scale matrix S = nu * sigma2
  m0 <- 12.5
  sigma0_sq <- 9.0
  kappa0 <- 4.0
  nu0 <- 8.0
  S0 <- matrix(nu0 * sigma0_sq, nrow = 1, ncol = 1, dimnames = list("VOT", "VOT"))

  nix_rep <- new_nix_category_representation(
    category_labels = "/b/", cue_labels = "VOT",
    m = m0, sigma2 = sigma0_sq, kappa = kappa0, nu = nu0
  )
  niw_rep <- new_niw_category_representation(
    category_labels = "/b/", cue_labels = "VOT",
    m = c(VOT = m0), S = S0, kappa = kappa0, nu = nu0
  )

  # 1. Parity for update_category_representation
  obs_N <- 5
  obs_mean <- c(VOT = 22.0)
  obs_css <- matrix(45.0, nrow = 1, ncol = 1, dimnames = list("VOT", "VOT"))

  up_nix_rep <- update_category_representation(nix_rep, x_N = obs_N, x_mean = obs_mean, x_css = obs_css)
  up_niw_rep <- update_category_representation(niw_rep, x_N = obs_N, x_mean = obs_mean, x_css = obs_css)

  expect_equal(unname(up_nix_rep@m), unname(up_niw_rep@m))
  expect_equal(up_nix_rep@kappa, up_niw_rep@kappa)
  expect_equal(up_nix_rep@nu, up_niw_rep@nu)
  # Mathematical identity: nu_n * sigma2_n == S_n
  expect_equal(up_nix_rep@nu * up_nix_rep@sigma2, unname(up_niw_rep@S[1, 1]))

  # 2. Parity for full models (update_template) across multiple categories
  cats <- c("/b/", "/d/")
  m_vec <- c(10, 45)
  sig2_vec <- c(4, 16)
  kap_vec <- c(3, 5)
  nu_vec <- c(10, 12)

  nix_reps <- list(
    new_nix_category_representation("/b/", "VOT", m = m_vec[1], sigma2 = sig2_vec[1], kappa = kap_vec[1], nu = nu_vec[1]),
    new_nix_category_representation("/d/", "VOT", m = m_vec[2], sigma2 = sig2_vec[2], kappa = kap_vec[2], nu = nu_vec[2])
  )
  names(nix_reps) <- cats
  nix_mod <- new_nix_ideal_adaptor(new_category_representation_template(nix_reps), category_prior = c(0.5, 0.5))

  niw_reps <- list(
    new_niw_category_representation("/b/", "VOT", m = c(VOT = m_vec[1]), S = matrix(nu_vec[1] * sig2_vec[1], dimnames = list("VOT", "VOT")), kappa = kap_vec[1], nu = nu_vec[1]),
    new_niw_category_representation("/d/", "VOT", m = c(VOT = m_vec[2]), S = matrix(nu_vec[2] * sig2_vec[2], dimnames = list("VOT", "VOT")), kappa = kap_vec[2], nu = nu_vec[2])
  )
  names(niw_reps) <- cats
  niw_mod <- new_niw_ideal_adaptor(new_category_representation_template(niw_reps), category_prior = c(0.5, 0.5))

  exposure <- tibble::tibble(
    category = c("/b/", "/b/", "/d/", "/b/", "/d/", "/d/"),
    VOT = c(8, 14, 40, 11, 48, 52)
  )

  # Batch update parity
  up_nix_batch <- update_template(nix_mod, exposure, updating = "batch")
  up_niw_batch <- update_template(niw_mod, exposure, updating = "batch")

  for (cat_name in cats) {
    r_nix <- up_nix_batch@category_template@representations[[cat_name]]
    r_niw <- up_niw_batch@category_template@representations[[cat_name]]
    expect_equal(unname(r_nix@m), unname(r_niw@m))
    expect_equal(r_nix@kappa, r_niw@kappa)
    expect_equal(r_nix@nu, r_niw@nu)
    expect_equal(r_nix@nu * r_nix@sigma2, unname(r_niw@S[1, 1]))
  }

  # Incremental update parity
  up_nix_inc <- update_template(nix_mod, exposure, updating = "incremental")
  up_niw_inc <- update_template(niw_mod, exposure, updating = "incremental")

  for (cat_name in cats) {
    r_nix_inc <- up_nix_inc@category_template@representations[[cat_name]]
    r_niw_inc <- up_niw_inc@category_template@representations[[cat_name]]
    expect_equal(unname(r_nix_inc@m), unname(r_niw_inc@m))
    expect_equal(r_nix_inc@kappa, r_niw_inc@kappa)
    expect_equal(r_nix_inc@nu, r_niw_inc@nu)
    expect_equal(r_nix_inc@nu * r_nix_inc@sigma2, unname(r_niw_inc@S[1, 1]))
  }
})

test_that("update_category_representation preserves representation when x_N == 0", {
  nix_rep <- example_nix_category_representation()
  expect_identical(
    update_category_representation(nix_rep, x_N = 0, x_mean = 20, x_css = 5),
    nix_rep
  )

  mnix_rep <- example_mnix_category_representation()
  expect_identical(
    update_category_representation(mnix_rep, x_N = 0, x_mean = c(10, 0), x_css = c(1, 1)),
    mnix_rep
  )

  niw_rep <- example_niw_category_representation(n_cues = 2)
  expect_identical(
    update_category_representation(niw_rep, x_N = 0, x_mean = c(10, 0), x_css = diag(2)),
    niw_rep
  )
})

test_that("update_category_representation errors informatively when x_N > 0 and inputs contain NA", {
  nix_rep <- example_nix_category_representation()
  expect_error(
    update_category_representation(nix_rep, x_N = 2, x_mean = NA_real_, x_css = 5),
    "x_mean cannot contain NA values when x_N > 0"
  )
  expect_error(
    update_category_representation(nix_rep, x_N = 2, x_mean = 10, x_css = NA_real_),
    "x_css cannot contain NA values when x_N > 0"
  )

  mnix_rep <- example_mnix_category_representation()
  expect_error(
    update_category_representation(mnix_rep, x_N = 1, x_mean = c(10, NA_real_), x_css = c(1, 1)),
    "x_mean cannot contain NA values when x_N > 0"
  )
  expect_error(
    update_category_representation(mnix_rep, x_N = 1, x_mean = c(10, 20), x_css = c(1, NA_real_)),
    "x_css cannot contain NA values when x_N > 0"
  )

  niw_rep <- example_niw_category_representation(n_cues = 2)
  expect_error(
    update_category_representation(niw_rep, x_N = 1, x_mean = c(NA_real_, 20), x_css = diag(2)),
    "x_mean cannot contain NA values when x_N > 0"
  )
  expect_error(
    update_category_representation(niw_rep, x_N = 1, x_mean = c(10, 20), x_css = matrix(c(1, NA_real_, 0, 1), 2, 2)),
    "x_css cannot contain NA values when x_N > 0"
  )
})


