
.io <- suppressMessages(suppressWarnings(example_mvg_ideal_observer(n_cues = 3)))
.cues <- get_cue_labels(.io)
.data <- sample_observations(.io, n = 50)

test_that("transform_cues - input (data)", {
  expect_error(
    transform_cues(
      cues = .cues[1],
      return.transformed.data = F,
      return.transform.parameters = F,
      return.transform.function = F,
      return.untransform.function = F))
  expect_error(
    transform_cues(
      data = NA,
      cues = .cues[1],
      return.transformed.data = F,
      return.transform.parameters = F,
      return.transform.function = F,
      return.untransform.function = F))
  expect_error(
    transform_cues(
      data = NULL,
      cues = .cues[1],
      return.transformed.data = F,
      return.transform.parameters = F,
      return.transform.function = F,
      return.untransform.function = F))
  expect_error(
    transform_cues(
      data = "data",
      cues = .cues[1],
      return.transformed.data = F,
      return.transform.parameters = F,
      return.transform.function = F,
      return.untransform.function = F))
})

test_that("transform_cues - input (cues)", {
  expect_no_error(
    transform_cues(
      data = .data,
      cues = .cues[1],
      return.transformed.data = F,
      return.transform.parameters = F,
      return.transform.function = F,
      return.untransform.function = F))
  expect_no_error(
    transform_cues(
      data = .data,
      cues = .cues[1],
      return.transformed.data = T,
      return.transform.parameters = F,
      return.transform.function = F,
      return.untransform.function = F))
  expect_no_error(
    transform_cues(
      data = .data,
      cues = .cues[1:2],
      return.transformed.data = F,
      return.transform.parameters = T,
      return.transform.function = F,
      return.untransform.function = F))
  expect_error(
    transform_cues(
      data = .data,
      return.transformed.data = F,
      return.transform.parameters = F,
      return.transform.function = T,
      return.untransform.function = F))
  expect_error(
    transform_cues(
      data = .data,
      cues = cue1,
      return.transformed.data = F,
      return.transform.parameters = F,
      return.transform.function = T,
      return.untransform.function = F))
  expect_error(
    transform_cues(
      data = .data,
      cues = NA,
      return.transformed.data = F,
      return.transform.parameters = F,
      return.transform.function = T,
      return.untransform.function = F))
  expect_error(
    transform_cues(
      data = .data,
      cues = NULL,
      return.transformed.data = F,
      return.transform.parameters = F,
      return.transform.function = T,
      return.untransform.function = F))
  expect_error(
    transform_cues(
      data = .data,
      cues = c("not a cue"),
      return.transformed.data = F,
      return.transform.parameters = F,
      return.transform.function = T,
      return.untransform.function = F))
})

test_that("transform_cues - input (center, scale, pca, attach)", {
  expect_no_error(
    transform_cues(
      data = .data,
      cues = .cues[1],
      center = F,
      scale = F,
      pca = F,
      attach = F,
      return.transformed.data = F,
      return.transform.parameters = F,
      return.transform.function = F,
      return.untransform.function = F))
  expect_no_error(
    transform_cues(
      data = .data,
      cues = .cues[1],
      center = T,
      scale = F,
      pca = F,
      attach = F,
      return.transformed.data = F,
      return.transform.parameters = F,
      return.transform.function = F,
      return.untransform.function = F))
  expect_no_error(
    transform_cues(
      data = .data,
      cues = .cues[1],
      center = F,
      scale = T,
      pca = F,
      attach = F,
      return.transformed.data = F,
      return.transform.parameters = F,
      return.transform.function = F,
      return.untransform.function = F))
  expect_no_error(
    transform_cues(
      data = .data,
      cues = .cues[1],
      center = F,
      scale = F,
      pca = T,
      attach = F,
      return.transformed.data = F,
      return.transform.parameters = F,
      return.transform.function = F,
      return.untransform.function = F))
  expect_no_error(
    transform_cues(
      data = .data,
      cues = .cues[1],
      center = F,
      scale = F,
      pca = F,
      attach = T,
      return.transformed.data = F,
      return.transform.parameters = F,
      return.transform.function = F,
      return.untransform.function = F))
  expect_no_error(
    transform_cues(
      data = .data,
      cues = .cues[1],
      center = T,
      scale = T,
      pca = F,
      attach = F,
      return.transformed.data = F,
      return.transform.parameters = F,
      return.transform.function = F,
      return.untransform.function = F))
  expect_no_error(
    transform_cues(
      data = .data,
      cues = .cues[1],
      center = T,
      scale = T,
      pca = T,
      attach = F,
      return.transformed.data = F,
      return.transform.parameters = F,
      return.transform.function = F,
      return.untransform.function = F))
  expect_no_error(
    transform_cues(
      data = .data,
      cues = .cues[1],
      center = T,
      scale = T,
      pca = T,
      attach = T,
      return.transformed.data = F,
      return.transform.parameters = F,
      return.transform.function = F,
      return.untransform.function = F))
  expect_error(
    transform_cues(
      data = .data,
      cues = .cues[1],
      center = "yes",
      scale = F,
      pca = F,
      attach = F,
      return.transformed.data = F,
      return.transform.parameters = F,
      return.transform.function = F,
      return.untransform.function = F))
  expect_error(
    transform_cues(
      data = .data,
      cues = .cues[1],
      center = F,
      scale = "yes",
      pca = F,
      attach = F,
      return.transformed.data = F,
      return.transform.parameters = F,
      return.transform.function = F,
      return.untransform.function = F))
  expect_error(
    transform_cues(
      data = .data,
      cues = .cues[1],
      center = F,
      scale = F,
      pca = "yes",
      attach = F,
      return.transformed.data = F,
      return.transform.parameters = F,
      return.transform.function = F,
      return.untransform.function = F))
  expect_error(
    transform_cues(
      data = .data,
      cues = .cues[1],
      center = F,
      scale = F,
      pca = F,
      attach = "yes",
      return.transformed.data = F,
      return.transform.parameters = F,
      return.transform.function = F,
      return.untransform.function = F))
})

test_that("transform_cues - input (return arguments)", {
  expect_no_error(
    transform_cues(
      data = .data,
      cues = .cues[1],
      return.transformed.data = F,
      return.transform.parameters = F,
      return.transform.function = F,
      return.untransform.function = F))
  expect_no_error(
    transform_cues(
      data = .data,
      cues = .cues[1],
      return.transformed.data = T,
      return.transform.parameters = F,
      return.transform.function = F,
      return.untransform.function = F))
  expect_no_error(
    transform_cues(
      data = .data,
      cues = .cues[1],
      return.transformed.data = F,
      return.transform.parameters = T,
      return.transform.function = F,
      return.untransform.function = F))
  expect_no_error(
    transform_cues(
      data = .data,
      cues = .cues[1],
      return.transformed.data = F,
      return.transform.parameters = F,
      return.transform.function = T,
      return.untransform.function = F))
  expect_no_error(
    transform_cues(
      data = .data,
      cues = .cues[1],
      return.transformed.data = F,
      return.transform.parameters = F,
      return.transform.function = F,
      return.untransform.function = T))
  expect_no_error(
    transform_cues(
      data = .data,
      cues = .cues[1],
      return.transformed.data = T,
      return.transform.parameters = T,
      return.transform.function = F,
      return.untransform.function = F))
  expect_no_error(
    transform_cues(
      data = .data,
      cues = .cues[1],
      return.transformed.data = T,
      return.transform.parameters = F,
      return.transform.function = T,
      return.untransform.function = F))
  expect_no_error(
    transform_cues(
      data = .data,
      cues = .cues[1],
      return.transformed.data = T,
      return.transform.parameters = F,
      return.transform.function = F,
      return.untransform.function = T))
  expect_no_error(
    transform_cues(
      data = .data,
      cues = .cues[1],
      return.transformed.data = T,
      return.transform.parameters = T,
      return.transform.function = T,
      return.untransform.function = F))
  expect_no_error(
    transform_cues(
      data = .data,
      cues = .cues[1],
      return.transformed.data = T,
      return.transform.parameters = F,
      return.transform.function = T,
      return.untransform.function = T))
  expect_no_error(
    transform_cues(
      data = .data,
      cues = .cues[1],
      return.transformed.data = T,
      return.transform.parameters = T,
      return.transform.function = T,
      return.untransform.function = T))
  expect_no_error(
    transform_cues(
      data = .data,
      cues = .cues[1],
      return.transformed.data = F,
      return.transform.parameters = T,
      return.transform.function = T,
      return.untransform.function = F))
  expect_no_error(
    transform_cues(
      data = .data,
      cues = .cues[1],
      return.transformed.data = F,
      return.transform.parameters = F,
      return.transform.function = T,
      return.untransform.function = T))
  expect_no_error(
    transform_cues(
      data = .data,
      cues = .cues[1],
      return.transformed.data = F,
      return.transform.parameters = T,
      return.transform.function = T,
      return.untransform.function = T))
})

test_that("transform_cues - basic output", {
  expect_message(
    transform_cues(
      data = .data,
      cues = .cues[1],
      return.transformed.data = F,
      return.transform.parameters = F,
      return.transform.function = F,
      return.untransform.function = F))
  expect_true(
    is.data.frame(
      transform_cues(
        data = .data,
        cues = .cues[1],
        return.transformed.data = T,
        return.transform.parameters = F,
        return.transform.function = F,
        return.untransform.function = F)))
  expect_true(
    is.list(
      transform_cues(
        data = .data,
        cues = .cues[1],
        center = T,
        scale = T,
        pca = F,
        return.transformed.data = F,
        return.transform.parameters = T,
        return.transform.function = F,
        return.untransform.function = F)))
  expect_true(
    is.function(
    transform_cues(
      data = .data,
      cues = .cues[1],
      center = T,
      scale = T,
      pca = F,
      return.transformed.data = F,
      return.transform.parameters = F,
      return.transform.function = T,
      return.untransform.function = F)))
  expect_true(
    is.function(
    transform_cues(
      data = .data,
      cues = .cues[1],
      center = T,
      scale = T,
      pca = F,
      return.transformed.data = F,
      return.transform.parameters = F,
      return.transform.function = F,
      return.untransform.function = T)))
  expect_true(
    is.list(
      transform_cues(
        data = .data,
        cues = .cues[1],
        center = T,
        scale = T,
        pca = F,
        return.transformed.data = T,
        return.transform.parameters = T,
        return.transform.function = F,
        return.untransform.function = F)))
  expect_true(
    is.list(
      transform_cues(
        data = .data,
        cues = .cues[1],
        center = T,
        scale = T,
        pca = F,
        return.transformed.data = T,
        return.transform.parameters = T,
        return.transform.function = T,
        return.untransform.function = F)))
  expect_true(
    is.list(
      transform_cues(
        data = .data,
        cues = .cues[1],
        center = T,
        scale = T,
        pca = F,
        return.transformed.data = T,
        return.transform.parameters = T,
        return.transform.function = T,
        return.untransform.function = T)))
  expect_true(
    is.data.frame(
      transform_cues(
        data = .data,
        cues = .cues[1],
        center = T,
        scale = T,
        pca = F,
        return.transformed.data = T,
        return.transform.parameters = T,
        return.transform.function = T,
        return.untransform.function = T)[["data"]]))
  expect_true(
    is.list(
      transform_cues(
        data = .data,
        cues = .cues[1],
        center = T,
        scale = T,
        pca = F,
        return.transformed.data = T,
        return.transform.parameters = T,
        return.transform.function = T,
        return.untransform.function = T)[["transform.parameters"]]))
  expect_true(
    is.function(
    transform_cues(
      data = .data,
      cues = .cues[1],
      center = T,
      scale = T,
      pca = F,
      return.transformed.data = T,
      return.transform.parameters = T,
      return.transform.function = T,
      return.untransform.function = T)[["transform.function"]]))
  expect_true(
    is.function(
    transform_cues(
      data = .data,
      cues = .cues[1],
      center = T,
      scale = T,
      pca = F,
      return.transformed.data = T,
      return.transform.parameters = T,
      return.transform.function = T,
      return.untransform.function = T)[["untransform.function"]]))
})

.data.transformed <- transform_cues(data = .data, cues = .cues, center = T, scale = T)
.transform <- transform_cues(data = .data, cues = .cues, center = T, scale = T, return.transformed.data = F, return.transform.parameters = T)
test_that("transform_cues - do transformed cue values have changed?", {
  expect_false(all(.data[[.cues[1]]] == .data.transformed[[.cues[1]]]))
  expect_false(all(.data[[.cues[2]]] == .data.transformed[[.cues[2]]]))
  expect_false(all(.data[[.cues[3]]] == .data.transformed[[.cues[3]]]))
})

test_that("transform_cues - output equality of input and un-transformation", {
  # one cue, without PCA
  tc <-
    transform_cues(
      data = .data,
      cues = .cues[1],
      center = T,
      scale = T,
      pca = F,
      return.transformed.data = F,
      return.transform.parameters = F,
      return.transform.function = T,
      return.untransform.function = T)
  expect_equal(tc[["untransform.function"]](tc[["transform.function"]](.data))[[.cues[1]]], .data[[.cues[1]]])
  tc <-
    transform_cues(
      data = .data,
      cues = .cues[1],
      center = T,
      scale = T,
      pca = F,
      return.transformed.data = T,
      return.transform.parameters = F,
      return.transform.function = T,
      return.untransform.function = T)
  expect_equal(tc[["untransform.function"]](tc[["transform.function"]](.data))[[.cues[1]]], .data[[.cues[1]]])
  tc <-
    transform_cues(
      data = .data,
      cues = .cues[1],
      center = T,
      scale = T,
      pca = F,
      return.transformed.data = F,
      return.transform.parameters = T,
      return.transform.function = T,
      return.untransform.function = T)
  expect_equal(tc[["untransform.function"]](tc[["transform.function"]](.data))[[.cues[1]]], .data[[.cues[1]]])
  tc <-
    transform_cues(
      data = .data,
      cues = .cues[1],
      center = T,
      scale = T,
      pca = F,
      return.transformed.data = T,
      return.transform.parameters = T,
      return.transform.function = T,
      return.untransform.function = T)
  expect_equal(tc[["untransform.function"]](tc[["transform.function"]](.data))[[.cues[1]]], .data[[.cues[1]]])
  # two cues, with PCA
  tc <-
    transform_cues(
      data = .data,
      cues = .cues[1:2],
      center = T,
      scale = T,
      pca = T,
      return.transformed.data = F,
      return.transform.parameters = F,
      return.transform.function = T,
      return.untransform.function = T)
  expect_equal(tc[["untransform.function"]](tc[["transform.function"]](.data))[[.cues[1]]], .data[[.cues[1]]])
  expect_equal(tc[["untransform.function"]](tc[["transform.function"]](.data))[[.cues[2]]], .data[[.cues[2]]])
  tc <-
    transform_cues(
      data = .data,
      cues = .cues[1:2],
      center = T,
      scale = T,
      pca = T,
      return.transformed.data = F,
      return.transform.parameters = T,
      return.transform.function = T,
      return.untransform.function = T)
  expect_equal(tc[["untransform.function"]](tc[["transform.function"]](.data))[[.cues[1]]], .data[[.cues[1]]])
  expect_equal(tc[["untransform.function"]](tc[["transform.function"]](.data))[[.cues[2]]], .data[[.cues[2]]])
  tc <-
    transform_cues(
      data = .data,
      cues = .cues[1:2],
      center = T,
      scale = T,
      pca = T,
      return.transformed.data = T,
      return.transform.parameters = F,
      return.transform.function = T,
      return.untransform.function = T)
  expect_equal(tc[["untransform.function"]](tc[["transform.function"]](.data))[[.cues[1]]], .data[[.cues[1]]])
  expect_equal(tc[["untransform.function"]](tc[["transform.function"]](.data))[[.cues[2]]], .data[[.cues[2]]])
  tc <-
    transform_cues(
      data = .data,
      cues = .cues[1:2],
      center = T,
      scale = T,
      pca = T,
      return.transformed.data = T,
      return.transform.parameters = T,
      return.transform.function = T,
      return.untransform.function = T)
  expect_equal(tc[["untransform.function"]](tc[["transform.function"]](.data))[[.cues[1]]], .data[[.cues[1]]])
  expect_equal(tc[["untransform.function"]](tc[["transform.function"]](.data))[[.cues[2]]], .data[[.cues[2]]])
})

test_that("untransform_mean back-transforms category means accurately", {
  for (n_cue in 1:3) {
    cues_sub <- .cues[seq_len(n_cue)]
    tr <- get_affine_transform(.data, cues = cues_sub, type = "standardize")
    trans_data <- tr@transform.function(.data, return_type = "replace")

    orig_means <- .data %>%
      get_sufficient_category_statistics(cues = cues_sub) %>%
      pull(mean)

    trans_means <- trans_data %>%
      get_sufficient_category_statistics(cues = cues_sub) %>%
      pull(mean)

    recovered_means <- purrr::map(trans_means, ~ untransform_mean(.x, tr))
    expect_equal(recovered_means, orig_means, tolerance = 1e-6)
  }
})

test_that("untransform_cov back-transforms category covariances accurately", {
  for (n_cue in 1:3) {
    cues_sub <- .cues[seq_len(n_cue)]
    tr <- get_affine_transform(.data, cues = cues_sub, type = "standardize")
    trans_data <- tr@transform.function(.data, return_type = "replace")

    orig_covs <- .data %>%
      get_sufficient_category_statistics(cues = cues_sub) %>%
      pull(cov)

    trans_covs <- trans_data %>%
      get_sufficient_category_statistics(cues = cues_sub) %>%
      pull(cov)

    recovered_covs <- purrr::map(trans_covs, ~ untransform_cov(.x, tr))
    expect_equal(recovered_covs, orig_covs, tolerance = 1e-6)
  }
})

# =============================================================================
# S7 Affine Transformations & Algebraic Inversion
# =============================================================================

test_that("get_affine_transform creates valid MVBU_TransformInformation", {
  # 1D transform
  tr1 <- get_affine_transform(.data, cues = .cues[1], type = "standardize")
  expect_true(S7::S7_inherits(tr1, MVBU_TransformInformation))
  expect_equal(tr1@cue.labels, .cues[1])
  expect_equal(tr1@type, "standardize")
  expect_length(tr1@shift, 1L)
  expect_true(is.numeric(tr1@SCALE))
  expect_false(is.null(tr1@transform.function))

  # Multi-D transform
  tr3 <- get_affine_transform(.data, cues = .cues, type = "PCA whiten")
  expect_true(S7::S7_inherits(tr3, MVBU_TransformInformation))
  expect_equal(tr3@cue.labels, .cues)
  expect_equal(tr3@type, "PCA whiten")
  expect_length(tr3@shift, 3L)
  expect_true(is.matrix(tr3@SCALE))
  expect_equal(dim(tr3@SCALE), c(3L, 3L))

  # Test transform.function on data
  d_trans <- tr3@transform.function(.data, return_type = "replace")
  cov_trans <- cov(as.matrix(d_trans[, .cues]))
  expect_equal(unname(diag(cov_trans)), rep(1, 3L), tolerance = 1e-4)
})

test_that("invert_transform satisfies algebraic duality and involution", {
  # 1D case
  tr1 <- get_affine_transform(.data, cues = .cues[1], type = "standardize")
  inv1 <- invert_transform(tr1)
  expect_true(S7::S7_inherits(inv1, MVBU_TransformInformation))
  expect_equal(inv1@SCALE, 1 / as.numeric(tr1@SCALE))
  expect_equal(inv1@shift, -as.numeric(tr1@SCALE) * tr1@shift)

  # Double inversion recovers original transform
  inv_inv1 <- invert_transform(inv1)
  expect_equal(inv_inv1@shift, tr1@shift, tolerance = 1e-7)
  expect_equal(inv_inv1@SCALE, tr1@SCALE, tolerance = 1e-7)

  # Multi-D case
  tr3 <- get_affine_transform(.data, cues = .cues, type = "ZCA whiten")
  inv3 <- invert_transform(tr3)
  expect_true(S7::S7_inherits(inv3, MVBU_TransformInformation))
  expect_equal(inv3@SCALE, solve(tr3@SCALE), tolerance = 1e-7)
  expect_equal(
    inv3@shift,
    -as.vector(tr3@SCALE %*% tr3@shift),
    tolerance = 1e-7
  )

  # Double inversion
  inv_inv3 <- invert_transform(inv3)
  expect_equal(inv_inv3@shift, tr3@shift, tolerance = 1e-7)
  expect_equal(inv_inv3@SCALE, tr3@SCALE, tolerance = 1e-7)
})

test_that("transform_mean and untransform_mean round-trip", {
  tr3 <- get_affine_transform(.data, cues = .cues, type = "standardize")
  m_orig <- c(cue1 = 10, cue2 = -5, cue3 = 2.5)

  m_trans <- transform_mean(m_orig, tr3)
  expect_named(m_trans, names(m_orig))
  expect_false(all(m_trans == m_orig))

  m_back <- untransform_mean(m_trans, tr3)
  expect_equal(m_back, m_orig, tolerance = 1e-7)

  # 1D round-trip
  tr1 <- get_affine_transform(.data, cues = .cues[1], type = "standardize")
  m1_orig <- c(cue1 = 4.2)
  m1_trans <- transform_mean(m1_orig, tr1)
  m1_back <- untransform_mean(m1_trans, tr1)
  expect_equal(m1_back, m1_orig, tolerance = 1e-7)
})

test_that("transform_cov and untransform_cov round-trip", {
  tr3 <- get_affine_transform(.data, cues = .cues, type = "PCA whiten")
  S_orig <- matrix(c(2, 0.5, 0.1, 0.5, 3, 0.2, 0.1, 0.2, 1.5), 3, 3)
  dimnames(S_orig) <- list(.cues, .cues)

  S_trans <- transform_cov(S_orig, tr3)
  expect_equal(dimnames(S_trans), dimnames(S_orig))

  S_back <- untransform_cov(S_trans, tr3)
  expect_equal(S_back, S_orig, tolerance = 1e-7)

  # 1D round-trip
  tr1 <- get_affine_transform(.data, cues = .cues[1], type = "standardize")
  S1_orig <- matrix(4, 1, 1, dimnames = list(.cues[1], .cues[1]))
  S1_trans <- transform_cov(S1_orig, tr1)
  S1_back <- untransform_cov(S1_trans, tr1)
  expect_equal(S1_back, S1_orig, tolerance = 1e-7)
})

test_that(
  "transform_model and untransform_model round-trip across representations", {
  tr <- get_affine_transform(.data, cues = .cues, type = "standardize")

  # MVG model
  mvg_model <- example_mvg_ideal_observer(n_cues = 3)
  mvg_trans <- transform_model(mvg_model, tr)
  mvg_back <- untransform_model(mvg_trans, tr)
  for (cat in names(mvg_model@category_template@representations)) {
    expect_equal(
      mvg_back@category_template@representations[[cat]]@mu,
      mvg_model@category_template@representations[[cat]]@mu,
      tolerance = 1e-6
    )
    expect_equal(
      mvg_back@category_template@representations[[cat]]@Sigma,
      mvg_model@category_template@representations[[cat]]@Sigma,
      tolerance = 1e-6
    )
  }

  # NIW model
  niw_model <- example_niw_ideal_adaptor(n_cues = 3)
  niw_trans <- transform_model(niw_model, tr)
  niw_back <- untransform_model(niw_trans, tr)
  for (cat in names(niw_model@category_template@representations)) {
    expect_equal(
      niw_back@category_template@representations[[cat]]@m,
      niw_model@category_template@representations[[cat]]@m,
      tolerance = 1e-6
    )
    expect_equal(
      niw_back@category_template@representations[[cat]]@S,
      niw_model@category_template@representations[[cat]]@S,
      tolerance = 1e-6
    )
  }

  # 1D UVG and NIX models
  tr1 <- get_affine_transform(.data, cues = .cues[1], type = "standardize")
  uvg_model <- example_uvg_ideal_observer()
  uvg_trans <- transform_model(uvg_model, tr1)
  uvg_back <- untransform_model(uvg_trans, tr1)
  for (cat in names(uvg_model@category_template@representations)) {
    expect_equal(
      uvg_back@category_template@representations[[cat]]@mu,
      uvg_model@category_template@representations[[cat]]@mu,
      tolerance = 1e-6
    )
    expect_equal(
      uvg_back@category_template@representations[[cat]]@sigma2,
      uvg_model@category_template@representations[[cat]]@sigma2,
      tolerance = 1e-6
    )
  }
})

test_that(
  "deprecated transform functions emit lifecycle deprecation warnings", {
  expect_warning(
    transform_cues(data = .data, cues = .cues[1]),
    class = "lifecycle_warning_deprecated"
  )
  expect_warning(
    untransform_cues(
      data = .data, cues = .cues[1],
      transform.parameters = list(cue.labels = .cues[1])
    ),
    class = "lifecycle_warning_deprecated"
  )
  tr1 <- get_affine_transform(.data, cues = .cues[1], type = "standardize")
  expect_warning(
    transform_category_mean(c(cue1 = 1), tr1),
    class = "lifecycle_warning_deprecated"
  )
  expect_warning(
    untransform_category_mean(c(cue1 = 1), tr1),
    class = "lifecycle_warning_deprecated"
  )
  expect_warning(
    transform_category_cov(matrix(1, 1, 1), tr1),
    class = "lifecycle_warning_deprecated"
  )
  expect_warning(
    untransform_category_cov(matrix(1, 1, 1), tr1),
    class = "lifecycle_warning_deprecated"
  )
})

