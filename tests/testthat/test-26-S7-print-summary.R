test_that("print and summary work on category representations", {
  mvg_rep <- new_mvg_category_representation(
    category_labels = "A",
    cue_labels = c("F1", "F2"),
    mu = c(0, 1),
    Sigma = diag(2)
  )

  out <- capture.output(print(mvg_rep))
  expect_true(any(grepl("MVG_CategoryRepresentation", out)))
  expect_true(any(grepl("Category: A", out)))
  expect_true(any(grepl("Cues \\(2\\): F1, F2", out)))

  out_sum <- capture.output(summary(mvg_rep))
  expect_equal(out, out_sum)
})

test_that("print on NIW displays kappa and nu before m and S", {
  niw_rep <- new_niw_category_representation(
    category_labels = "A",
    cue_labels = c("F1", "F2"),
    m = c(0, 1),
    S = diag(2),
    kappa = 10,
    nu = 5
  )
  out <- capture.output(print(niw_rep))
  kappa_idx <- grep("kappa:", out)
  nu_idx <- grep("nu:", out)
  m_idx <- grep("m:", out)
  S_idx <- grep("S:", out)
  expect_true(length(kappa_idx) == 1 && length(nu_idx) == 1 && length(m_idx) == 1 && length(S_idx) == 1)
  expect_true(kappa_idx < m_idx)
  expect_true(nu_idx < m_idx)
  expect_true(m_idx < S_idx)
})

test_that("print on Exemplar representation displays first 5 exemplars", {
  ex_mat <- matrix(seq(1, 20), ncol = 2)
  colnames(ex_mat) <- c("F1", "F2")
  ex_rep <- new_exemplar_category_representation(
    category_labels = "A",
    cue_labels = c("F1", "F2"),
    exemplars = ex_mat
  )
  out <- capture.output(print(ex_rep))
  expect_true(any(grepl("Exemplars \\(10 points\\):", out)))
  expect_true(any(grepl("\\[1\\]", out)))
  expect_true(any(grepl("\\[5\\]", out)))
  expect_true(any(grepl("\\.\\.\\.", out)))
})

test_that("print and summary work on category representation templates", {
  mvg_rep <- new_mvg_category_representation(
    category_labels = "A",
    cue_labels = c("F1", "F2"),
    mu = c(0, 1),
    Sigma = diag(2)
  )
  tpl <- new_category_representation_template(representations = list(A = mvg_rep))

  out <- capture.output(print(tpl))
  expect_true(any(grepl("CategoryRepresentationTemplate", out)))
  expect_true(any(grepl("Categories \\(1\\): A", out)))
  expect_true(any(grepl("\\$A: MVG\\(mu = vector\\(2\\), Sigma = matrix\\(2, 2\\)\\)", out)))

  out_sum <- capture.output(summary(tpl))
  expect_equal(out, out_sum)
})

test_that("print and summary work on cognitive models", {
  mvg_rep <- new_mvg_category_representation(
    category_labels = "A",
    cue_labels = c("F1", "F2"),
    mu = c(0, 1),
    Sigma = diag(2)
  )
  tpl <- new_category_representation_template(representations = list(A = mvg_rep))
  model <- new_mvg_ideal_observer(
    category_template = tpl,
    category_prior = c(A = 1),
    lapse_rate = 0.05
  )

  out <- capture.output(print(model))
  expect_true(any(grepl("MVG_IdealObserver", out)))
  expect_true(any(grepl("Decision rule: sampling", out)))
  expect_true(any(grepl("Lapse rate: 0.05 \\(treatment: no_lapses\\)", out)))
  expect_true(any(grepl("Perceptual noise: none \\(treatment: no_noise\\)", out)))
  expect_true(any(grepl("\\$A: MVG\\(mu = vector\\(2\\), Sigma = matrix\\(2, 2\\)\\)", out)))

  out_sum <- capture.output(summary(model))
  expect_equal(out, out_sum)
})

test_that("summary works on IdealAdaptorStanfitInput", {
  stan_input <- get_example_staninput(1, stanmodel = "NIW_ideal_adaptor")
  sum_inp <- summary(stan_input)
  expect_true(S7::S7_inherits(sum_inp, Summary_IdealAdaptorStanfitInput))
  expect_true(is.data.frame(sum_inp@exposure_statistics))
  expect_true(sum_inp@test_summary$n_observations > 0)

  out <- capture.output(print(sum_inp))
  expect_true(any(grepl("Exposure sufficient statistics", out)))
  expect_true(any(grepl("Test data summary", out)))
})

test_that("print and summary work on MVBU_Stanfit", {
  fit_file <- testthat::test_path("models", "example-stanfit-NIW_ideal_adaptor-1-standardize-42.rds")
  skip_if_not(file.exists(fit_file), "Fixture file not found")
  fit <- readRDS(fit_file)
  out <- capture.output(print(fit))
  expect_true(any(grepl("MVBU_Stanfit", out) | grepl("IdealAdaptorStanfit", out)))
  expect_true(any(grepl("Model type: NIW", out)))
  expect_true(any(grepl("Categories", out)))
  expect_true(any(grepl("Draws:", out)))

  sum_obj <- summary(fit)
  expect_true(S7::S7_inherits(sum_obj, Summary_MVBU_Stanfit))
  sum_df <- as.data.frame(sum_obj)
  expect_true(is.data.frame(sum_df))
  expect_true("Parameter" %in% names(sum_df))
  expect_true("Group" %in% names(sum_df) || "Category" %in% names(sum_df))

  out_sum <- capture.output(print(sum_obj))
  expect_true(any(grepl("Fitted parameters", out_sum)))
})
