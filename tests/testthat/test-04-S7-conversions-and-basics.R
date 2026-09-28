test_that("to_array converts vectors into arrays with requested inner dimensions", {
  x <- to_array(c(1, 2), inner_dims = 2, outer_dims = 2, simplify = FALSE)

  expect_equal(dim(x), c(2, 2))
  expect_equal(x, matrix(c(1, 2, 1, 2), nrow = 2, byrow = TRUE))
})

test_that("to_array stacks list inputs into expected higher-dimensional array", {
  x <- to_array(
    list(
      matrix(1:4, nrow = 2, byrow = TRUE),
      matrix(5:8, nrow = 2, byrow = TRUE)
    ),
    inner_dims = c(2, 2),
    outer_dims = 2,
    simplify = FALSE
  )

  expect_equal(dim(x), c(2, 2, 2))
  expect_equal(x[, , 1], matrix(1:4, nrow = 2, byrow = TRUE))
  expect_equal(x[, , 2], matrix(5:8, nrow = 2, byrow = TRUE))
})

test_that("get_expected_mu_from_m and get_m_from_expected_mu work correctly", {
  m_val <- c(0, 1)
  expect_equal(get_expected_mu_from_m(m_val), m_val)
  expect_equal(get_m_from_expected_mu(m_val), m_val)

  m_list <- list(c(0, 1), c(2, 3))
  expect_equal(get_expected_mu_from_m(m_list), m_list)
  expect_equal(get_m_from_expected_mu(m_list), m_list)
})

test_that("get_expected_Sigma_from_S and get_S_from_expected_Sigma work correctly", {
  S_mat <- diag(c(4, 9))
  nu_val <- 10
  D <- 2
  # Sigma = S / (nu - D - 1) = S / 7
  expected_Sigma <- S_mat / 7
  expect_equal(get_expected_Sigma_from_S(S_mat, nu_val), expected_Sigma)
  expect_equal(get_S_from_expected_Sigma(expected_Sigma, nu_val), S_mat)
})

.io <- suppressMessages(suppressWarnings(example_mvg_ideal_observer(n_cues = 2)))
.cues <- get_cue_labels(.io)
.data <- sample_observations(.io, Ns = 50)

test_that("uss2css, css2cov - does sum-of-square to cov conversion work?", {
  stats <- get_sufficient_category_statistics(.data, cues = .cues)
  cov_direct <- stats$x_cov
  css_from_uss <- mapply(
    uss2css, stats$x_uss, stats$x_N, stats$x_mean,
    SIMPLIFY = FALSE
  )
  cov_from_uss <- mapply(css2cov, css_from_uss, stats$x_N, SIMPLIFY = FALSE)

  expect_equal(cov_direct, cov_from_uss, ignore_attr = TRUE)
})

test_that("css2uss and conversion functions handle edge cases", {
  mat <- matrix(
    c(4, 1, 1, 3), 2, 2,
    dimnames = list(c("x", "y"), c("x", "y"))
  )
  mu <- c(x = 1, y = 2)

  # N = 0 returns NA matrix with preserved dimensions and dimnames
  na_mat <- css2uss(mat, 0, mu)
  expect_true(is.matrix(na_mat))
  expect_equal(dim(na_mat), c(2, 2))
  expect_equal(dimnames(na_mat), dimnames(mat))
  expect_true(all(is.na(na_mat)))

  # N < 0 throws an error
  expect_error(css2uss(mat, -1, mu))
  expect_error(css2cov(mat, -1))
  expect_error(cov2css(mat, -1))
  expect_error(uss2css(mat, -1, mu))

  # Invertibility for N > 1
  uss_val <- css2uss(mat, 10, mu)
  css_rec <- uss2css(uss_val, 10, mu)
  expect_equal(css_rec, mat)

  # uss2cov returns cov
  cov_from_uss <- uss2cov(uss_val, 10, mu)
  expect_equal(cov_from_uss, css2cov(mat, 10))
})

test_that("get_sufficient_category_statistics handles categories and filtering", {
  df <- tibble::tibble(
    category = factor(c("A", "A", "A", "B", "B", "C")),
    c1 = c(1, 2, 3, 10, 20, 100),
    c2 = c(4, 5, 6, 40, 50, 400)
  )
  # Request categories c("A", "B")
  res <- get_sufficient_category_statistics(
    df,
    cues = c("c1", "c2"),
    category = "category",
    categories = c("A", "B")
  )
  expect_equal(nrow(res), 2)
  expect_equal(as.character(res$category), c("A", "B"))
  expect_equal(res$x_N, c(3, 2))

  # Check category A mean and covariance
  expect_equal(unname(res$x_mean[[1]]), c(2, 5))
  expect_equal(dim(res$x_css[[1]]), c(2, 2))
  expect_equal(dim(res$x_cov[[1]]), c(2, 2))
})

test_that("get_sufficient_category_statistics computes family statistics", {
  df <- tibble::tibble(
    category = factor(c("A", "A", "B", "B")),
    c1 = c(1, 3, 10, 20),
    c2 = c(2, 4, 30, 40)
  )

  # NIX scalar sums of squares
  res_nix <- get_sufficient_category_statistics(
    df,
    cues = "c1",
    model_family = "NIX"
  )
  expect_true(is.numeric(res_nix$x_ss[[1]]))
  expect_equal(length(res_nix$x_ss[[1]]), 1L)
  expect_false("x_cov" %in% names(res_nix))
  expect_equal(res_nix$x_ss[[1]], sum((c(1, 3) - 2)^2))

  # MNIX vector sums of squares
  res_mnix <- get_sufficient_category_statistics(
    df,
    cues = c("c1", "c2"),
    model_family = "MNIX"
  )
  expect_equal(length(res_mnix$x_ss[[1]]), 2L)
  expect_equal(names(res_mnix$x_ss[[1]]), c("c1", "c2"))
  expect_false("x_cov" %in% names(res_mnix))

  # NIW matrix sums of squares and covariance
  res_niw <- get_sufficient_category_statistics(
    df,
    cues = c("c1", "c2"),
    model_family = "NIW"
  )
  expect_true(is.matrix(res_niw$x_ss[[1]]))
  expect_true(is.matrix(res_niw$x_cov[[1]]))
  expect_equal(dim(res_niw$x_ss[[1]]), c(2, 2))
})

test_that("get_sufficient_category_statistics dispatches on MVBU_CognitiveModel", {
  io <- suppressMessages(suppressWarnings(
    example_mvg_ideal_observer(n_cues = 2)
  ))
  suff <- get_sufficient_category_statistics(io)
  expect_true(is.data.frame(suff))
  expect_true("category" %in% names(suff))
  expect_true("x_mean" %in% names(suff))
  expect_true("x_ss" %in% names(suff))
  expect_equal(nrow(suff), length(get_category_labels(io)))
})

test_that(".compute_default_cue_limits works for models, templates, and frames", {
  io <- suppressMessages(suppressWarnings(
    example_mvg_ideal_observer(n_cues = 2)
  ))
  cues <- get_cue_labels(io)

  lims_model <- .compute_default_cue_limits(io, cues = cues, n_sds = 2.0)
  expect_equal(names(lims_model), cues)
  expect_equal(length(lims_model[[1]]), 2L)
  expect_lt(lims_model[[1]][1], lims_model[[1]][2])

  lims_tmpl <- .compute_default_cue_limits(
    io@category_template,
    cues = cues,
    n_sds = 2.0
  )
  expect_equal(lims_model, lims_tmpl)
})

test_that("deprecated basics functions warn when called", {
  df <- data.frame(c1 = c(1, 2, 3), c2 = c(4, 5, 6))
  expect_warning(
    ss <- get_sum_of_squares_from_df(df, variables = c("c1", "c2")),
    class = "lifecycle_warning_deprecated"
  )
  expect_true(is.matrix(ss))

  expect_warning(
    vec_df <- make_vector_column(df, c("c1", "c2"), "vec"),
    class = "lifecycle_warning_deprecated"
  )
  expect_true("vec" %in% names(vec_df))
})
